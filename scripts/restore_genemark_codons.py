#!/usr/bin/env python3
"""
Put the start_codon and stop_codon rows of GeneMark genes back after
filterGenemark.pl.

filterGenemark.pl (add_single_cds) keeps only the CDS line of single-exon
genes, so in ET and EP mode every single-exon gene in genemark.f.good.gtf and
genemark.f.bad.gtf lost its start_codon and stop_codon rows. TSEBRA keeps
genemark.f.good.gtf unchanged (--keep_gtf), and these genes reached
braker.gtf without codon rows: no start codon in braker.gtf/braker.gff3, and
stringtie2utr labelled the stop codon as a 3 bp 3' UTR once normalize_gtf.py
had trimmed it from the CDS.

For each transcript in a filtered file that lacks start_codon or stop_codon
rows, the missing rows are copied from the GeneMark GTF that went into
filterGenemark.pl (same transcript_id, same sequence and strand, at an end of
the transcript's CDS). Rows are only ever added; nothing else changes.

Usage:
    restore_genemark_codons.py --source genemark.sorted.gtf \\
        genemark.f.good.gtf genemark.f.bad.gtf      (edited in place)
"""

import argparse
import os
import re
import sys
from collections import defaultdict

TX_RE = re.compile(r'transcript_id "([^"]+)"')
CODONS = ("start_codon", "stop_codon")


def tx_id(cols):
    m = TX_RE.search(cols[8]) if len(cols) >= 9 else None
    return m.group(1) if m else None


def read_codons(path):
    """transcript_id -> {codon type: [lines]} from the GeneMark GTF."""
    codons = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 9 and cols[2] in CODONS:
                tid = tx_id(cols)
                if tid:
                    codons[tid][cols[2]].append(line if line.endswith("\n") else line + "\n")
    return codons


def restore(path, codons):
    with open(path) as fh:
        lines = fh.readlines()
    order, by_tx, present = [], defaultdict(list), defaultdict(set)
    cds = defaultdict(list)
    for i, line in enumerate(lines):
        cols = line.rstrip("\n").split("\t")
        tid = tx_id(cols)
        if tid is None:
            continue
        if tid not in by_tx:
            order.append(tid)
        by_tx[tid].append(i)
        if cols[2] in CODONS:
            present[tid].add(cols[2])
        elif cols[2] == "CDS":
            cds[tid].append((cols[0], cols[6], int(cols[3]), int(cols[4])))

    insert_after = {}   # line index -> lines to add after it
    added = 0
    for tid in order:
        if tid not in cds:
            continue
        seqs = {(s, st) for s, st, _, _ in cds[tid]}
        if len(seqs) != 1:
            continue
        seqid, strand = seqs.pop()
        lo = min(a for _, _, a, _ in cds[tid])
        hi = max(b for _, _, _, b in cds[tid])
        new = []
        for ctype in CODONS:
            if ctype in present[tid]:
                continue
            for cl in codons.get(tid, {}).get(ctype, []):
                c = cl.rstrip("\n").split("\t")
                s, e = int(c[3]), int(c[4])
                # codon must belong to this CDS: same sequence and strand, at
                # one of its ends (GeneMark CDS include both codons)
                if c[0] != seqid or c[6] != strand or not (lo <= s and e <= hi):
                    continue
                if s != lo and e != hi:
                    continue
                new.append(cl)
        if new:
            insert_after[by_tx[tid][-1]] = new
            added += len(new)

    if added:
        out = []
        for i, line in enumerate(lines):
            out.append(line)
            out.extend(insert_after.get(i, []))
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            fh.writelines(out)
        os.replace(tmp, path)
    return added


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source", required=True,
                    help="GeneMark GTF given to filterGenemark.pl")
    ap.add_argument("gtf", nargs="+", help="filtered GTF files, edited in place")
    args = ap.parse_args()
    codons = read_codons(args.source)
    for path in args.gtf:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        n = restore(path, codons)
        print(f"[INFO] {os.path.basename(path)}: {n} start/stop codon rows restored "
              f"from {os.path.basename(args.source)}", file=sys.stderr)


if __name__ == "__main__":
    main()
