#!/usr/bin/env python3
"""
Drop-in replacement for ProtHint's bin/run_spliced_alignment.pl (issue #98).

The original dispatcher runs Spaln batches from Perl ithreads that share a
Thread::Queue and an MCE::Mutex and start each batch with system(). Users
report runs that stop at "Enqueueing pair ... (99.9%)" and sit idle for hours,
with braker.pl as well as BRAKER4. That progress line can never reach 100%
(the counter is printed before it is incremented), so the log looks the same
wherever the run stops once enqueueing is done. Enqueueing is throttled to
2 x cores batches ahead of the workers, so a single worker that gets stuck at
any point is invisible until the final join, which then waits forever.

This version keeps the command line, the temporary file names, the call to
spalnBatch.sh and the final gff_from_region_to_contig.pl step, but
  - runs each batch as a separate process group from a single-threaded
    process (no ithreads, no shared mutex, no fork from a threaded process),
  - kills a batch whose processes use no CPU time for --stall_timeout
    seconds, retries it once and then skips it with a warning,
  - writes the batch outputs in batch order, so spaln.gff does not depend on
    thread scheduling,
  - reports progress up to 100% and exits non-zero if spalnBatch.sh or the
    coordinate translation cannot run at all.

ProtHint finds the helper scripts next to this file: the directory of
sys.argv[0] without resolving symlinks (override with --prothint_bin or
PROTHINT_BIN). scripts/make_prothint_shadow.sh installs it into a symlinked
copy of the container's ETP bin tree.

Usage (as called by prothint.py):
    run_spliced_alignment.pl --cores N --nuc ../nuc.fasta --list pairs \
        --prot proteins.fa --v --aligner spaln --min_exon_score 25 \
        --longGene 30000 --longProtein 15000
"""

import argparse
import os
import signal
import subprocess
import sys
import time

BATCH_SIZE = 100
SPALN_OUT = "spaln.gff"
REGIONS_OUT = "spaln.regions.gff"
POLL_SECONDS = 0.2
CPU_CHECK_SECONDS = 30
KILL_WAIT_SECONDS = 30
MAX_ATTEMPTS = 2
DEFAULT_STALL_TIMEOUT = 1800

VERBOSE = False


def log(msg, force=False):
    if VERBOSE or force:
        sys.stderr.write("[" + time.ctime() + "] " + msg + "\n")
        sys.stderr.flush()


def read_list(path):
    """Read the nucleotide/protein pair list, same rules as ReadList."""
    pairs = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                sys.stderr.write("error, unexpected file format found "
                                 + sys.argv[0] + ": " + line)
                sys.exit(1)
            pairs.append((fields[0], fields[1]))
    return pairs


def read_sequences(path, wanted):
    """Read the FASTA records in `wanted`, uppercased and letters only,
    same rules as ReadSequence."""
    keep = bytes(range(ord("A"), ord("Z") + 1))
    drop = bytes(c for c in range(256) if c not in keep)
    seqs = {}
    current = None
    chunks = None
    have_id = False
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if current is not None:
                    seqs[current] = b"".join(chunks)
                fields = line[1:].split()
                name = fields[0].decode() if fields else ""
                have_id = bool(name)
                current = name if name in wanted else None
                chunks = []
            else:
                if not have_id:
                    sys.stderr.write("error, fasta record whithout definition "
                                     "line found " + sys.argv[0] + "\n")
                    sys.exit(1)
                if current is not None:
                    chunks.append(line.upper().translate(None, drop))
    if current is not None:
        seqs[current] = b"".join(chunks)
    return seqs


def save_single_fasta(path, seq_id, seq):
    with open(path, "wb") as out:
        out.write(b">" + seq_id.encode() + b"\n" + seq + b"\n")


def cpu_ticks_by_group():
    """CPU ticks used per process group, including reaped children.
    None if /proc cannot be read."""
    ticks = {}
    try:
        entries = os.listdir("/proc")
    except OSError:
        return None
    for entry in entries:
        if not entry.isdigit():
            continue
        try:
            with open("/proc/" + entry + "/stat") as fh:
                stat = fh.read()
        except OSError:
            continue
        # Fields after the command name, which may contain spaces:
        # [2] pgrp, [11:15] utime stime cutime cstime
        fields = stat[stat.rfind(")") + 2:].split()
        if len(fields) > 14:
            pgid = int(fields[2])
            used = sum(int(x) for x in fields[11:15])
            ticks[pgid] = ticks.get(pgid, 0) + used
    return ticks


class Batch:
    def __init__(self, index, pairs, first_counter):
        self.index = index
        self.pairs = pairs
        self.first_counter = first_counter
        self.name = "batch_" + str(index)
        self.out = self.name + "_out"
        self.attempts = 0
        self.proc = None
        self.cpu = None
        self.last_progress = 0.0

    def temp_files(self):
        for k in range(len(self.pairs)):
            counter = self.first_counter + k
            nuc = "nuc_" + str(counter)
            prot = "prot_" + str(counter)
            yield nuc, prot

    def write_inputs(self, nuc_seqs, prot_seqs):
        with open(self.name, "w") as batch_fh:
            for (nuc_id, prot_id), (nuc, prot) in zip(self.pairs,
                                                      self.temp_files()):
                save_single_fasta(nuc, nuc_id, nuc_seqs.get(nuc_id, b""))
                save_single_fasta(prot, prot_id, prot_seqs.get(prot_id, b""))
                batch_fh.write(nuc + "\t" + prot + "\n")

    def remove_leftovers(self):
        paths = [self.name, self.out]
        for nuc, prot in self.temp_files():
            paths += [nuc, prot, nuc + "_" + prot]
        for path in paths:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass


def kill_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    # A process in uninterruptible I/O only dies once the I/O returns;
    # do not let it block the dispatcher.
    try:
        proc.wait(timeout=KILL_WAIT_SECONDS)
    except subprocess.TimeoutExpired:
        log("WARNING: process group " + str(proc.pid) + " did not exit "
            "after SIGKILL, leaving it behind", force=True)


def append_file(out_fh, path):
    try:
        with open(path, "rb") as in_fh:
            while True:
                block = in_fh.read(1 << 20)
                if not block:
                    break
                out_fh.write(block)
    except FileNotFoundError:
        sys.stderr.write("Could not open '" + path + "' for reading\n")


def run(args):
    bin_dir = args.prothint_bin
    batch_script = os.path.join(bin_dir, "spalnBatch.sh")
    to_contig = os.path.join(bin_dir, "gff_from_region_to_contig.pl")
    for path in (batch_script, to_contig):
        if not os.access(path, os.X_OK):
            sys.stderr.write("error, ProtHint script not found or not "
                             "executable: " + path + "\n")
            sys.exit(1)

    os.environ["ALN_TAB"] = os.path.join(bin_dir, "..", "dependencies",
                                         "spaln_table")

    log("Starting spliced alignment with " + args.aligner)
    log("Loading alignment pairs into memory")
    pairs = read_list(args.list)
    nuc_seqs = read_sequences(args.nuc, {p[0] for p in pairs})
    prot_seqs = read_sequences(args.prot, {p[1] for p in pairs})
    if args.debug:
        for nuc_id, prot_id in pairs:
            if nuc_id not in nuc_seqs:
                sys.stderr.write(nuc_id + " missing\n")
            if prot_id not in prot_seqs:
                sys.stderr.write(prot_id + " missing\n")

    n_pairs = len(pairs)
    batches = [Batch(i // BATCH_SIZE + 1, pairs[i:i + BATCH_SIZE], i + 1)
               for i in range(0, n_pairs, BATCH_SIZE)]
    log("Pairs loaded. Number of pairs to align: " + str(n_pairs))
    log("Starting the alignments (" + str(len(batches)) + " batches of up to "
        + str(BATCH_SIZE) + " pairs, " + str(args.cores) + " in parallel, "
        + "stall timeout " + str(args.stall_timeout) + " s)")

    extra = ["{:g}".format(args.min_exon_score),
             "1" if args.nonCanonical else "0",
             str(args.longGene), str(args.longProtein)]

    pending = list(reversed(batches))
    running = {}
    finished = {}
    skipped = []
    next_to_write = 1
    n_done = 0
    next_print = 0
    start = last_cpu_check = time.time()

    def terminate(signum, frame):
        for batch in running.values():
            kill_group(batch.proc)
        sys.exit(128 + signum)

    signal.signal(signal.SIGTERM, terminate)
    signal.signal(signal.SIGINT, terminate)

    with open(REGIONS_OUT, "wb") as regions:
        while pending or running:
            while pending and len(running) < args.cores:
                batch = pending.pop()
                batch.attempts += 1
                batch.write_inputs(nuc_seqs, prot_seqs)
                try:
                    batch.proc = subprocess.Popen(
                        [batch_script, batch.name, batch.out] + extra,
                        start_new_session=True)
                except OSError as err:
                    sys.stderr.write("error, cannot run " + batch_script
                                     + ": " + str(err) + "\n")
                    for other in running.values():
                        kill_group(other.proc)
                    sys.exit(1)
                batch.cpu = None
                batch.last_progress = time.time()
                running[batch.index] = batch

            time.sleep(POLL_SECONDS)
            now = time.time()

            ticks = None
            if args.stall_timeout > 0 \
                    and now - last_cpu_check >= CPU_CHECK_SECONDS:
                last_cpu_check = now
                ticks = cpu_ticks_by_group()

            for index in list(running):
                batch = running[index]
                if batch.proc.poll() is not None:
                    del running[index]
                    finished[index] = batch
                    continue
                if ticks is None:
                    continue
                cpu = ticks.get(batch.proc.pid)
                if cpu is None or cpu != batch.cpu:
                    batch.cpu = cpu
                    batch.last_progress = now
                    continue
                if now - batch.last_progress < args.stall_timeout:
                    continue
                # No CPU time used by any process of this batch for
                # stall_timeout seconds: it is stuck, not slow.
                kill_group(batch.proc)
                del running[index]
                batch.remove_leftovers()
                if batch.attempts < MAX_ATTEMPTS:
                    log("WARNING: " + batch.name + " used no CPU for "
                        + str(args.stall_timeout) + " s, killed it and "
                        "retrying", force=True)
                    pending.append(batch)
                else:
                    log("WARNING: " + batch.name + " stalled again, skipping "
                        + "its " + str(len(batch.pairs)) + " pairs",
                        force=True)
                    skipped.append(batch)
                    finished[index] = None

            while next_to_write in finished:
                batch = finished.pop(next_to_write)
                if batch is not None:
                    append_file(regions, batch.out)
                    batch.remove_leftovers()
                    n_done += len(batch.pairs)
                next_to_write += 1
                regions.flush()

            permille = n_done * 1000 // n_pairs if n_pairs else 1000
            if permille >= next_print and n_done:
                elapsed = now - start
                seconds_left = int((1000 - permille) * elapsed / permille) \
                    if permille else 0
                log("Aligned pair %d/%d (%.1f%%). Est. time left: "
                    "%02d:%02d:%02d (hh:mm:ss)"
                    % (n_done, n_pairs, permille / 10, seconds_left // 3600,
                       seconds_left // 60 % 60, seconds_left % 60))
                next_print = permille + 1

    n_skipped = sum(len(b.pairs) for b in skipped)
    if n_skipped:
        log("WARNING: " + str(n_skipped) + " of " + str(n_pairs) + " pairs in "
            + str(len(skipped)) + " stalled batches were not aligned: "
            + ", ".join(b.name for b in skipped), force=True)
    log(str(n_pairs - n_skipped) + "/" + str(n_pairs) + " pairs aligned")
    log("Alignment of pairs finished")
    log("Translating coordinates from local pair level to contig level")

    rc = subprocess.call([to_contig, "--in_gff", REGIONS_OUT,
                          "--seq", args.nuc, "--out_gff", SPALN_OUT])
    if rc != 0:
        sys.stderr.write("error, " + to_contig + " exited with " + str(rc)
                         + "\n")
        sys.exit(1)
    os.unlink(REGIONS_OUT)
    log("Finished spliced alignment")


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Run Spaln on gene/protein pairs (ProtHint dispatcher "
                    "without Perl threads).")
    parser.add_argument("--nuc", required=True)
    parser.add_argument("--prot", required=True)
    parser.add_argument("--list", required=True)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--min_exon_score", type=float, default=25)
    parser.add_argument("--nonCanonical", action="store_true")
    parser.add_argument("--longGene", type=int, default=30000)
    parser.add_argument("--longProtein", type=int, default=15000)
    parser.add_argument(
        "--stall_timeout", type=int,
        default=int(os.environ.get("BRAKER4_SPALN_STALL_TIMEOUT",
                                   DEFAULT_STALL_TIMEOUT)),
        help="Kill a batch whose processes use no CPU for this many seconds "
             "(0 disables; default %(default)s, or "
             "BRAKER4_SPALN_STALL_TIMEOUT)")
    parser.add_argument(
        "--prothint_bin",
        default=os.environ.get(
            "PROTHINT_BIN", os.path.dirname(os.path.abspath(sys.argv[0]))),
        help="ProtHint bin directory with spalnBatch.sh (default: the "
             "directory this script was called from, symlinks unresolved)")
    args = parser.parse_args(argv)

    for name in ("nuc", "prot", "list"):
        path = getattr(args, name)
        if not os.path.exists(path):
            sys.stderr.write("error, file not found " + sys.argv[0] + ": "
                             + path + "\n")
            sys.exit(1)
        setattr(args, name, os.path.abspath(path))
    if args.cores < 1:
        sys.stderr.write("error, out of range number of cores: "
                         + str(args.cores) + "\n")
        sys.exit(1)
    args.aligner = args.aligner.lower()
    if args.aligner != "spaln":
        sys.stderr.write("error, invalid aligner specified: " + args.aligner
                         + ". Only \"Spaln\" is supported.\n")
        sys.exit(1)
    return args


def main(argv=None):
    global VERBOSE
    args = parse_args(sys.argv[1:] if argv is None else argv)
    VERBOSE = args.verbose or args.debug
    run(args)


if __name__ == "__main__":
    main()
