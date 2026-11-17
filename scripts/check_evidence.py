#!/usr/bin/env python3
"""Decide whether the evidence of a sample is enough for the planned
GeneMark mode (automatic mode fallback, [PARAMS] auto_mode_fallback).

Two subcommands, each writes a small decision file (key<TAB>value lines)
that rules/common.smk reads through get_effective_mode():

introns   Count the intron hints of a GFF file and those with a multiplicity
          of at least --min-mult. Fewer than --min-introns of the latter
          means: not enough for GeneMark-ET (RNA-Seq, multiplicity >= 10) or
          GeneMark-EP (proteins, multiplicity >= 4). The numbers are those
          of check_genemark_hints() in braker.pl, which warns that GeneMark
          may crash below 150 such introns.

result    Result of a GeneMark run (--tool GeneMark-ETP, GeneMark-ET or
          GeneMark-EP; 'etp' is an alias of this subcommand). The GeneMark
          rule writes a marker file when the program stopped for lack of
          training data (GeneMark-ETP: no spliced transcripts, no gene
          seeds, or too few training genes; GeneMark-ET/EP: an error of
          gmes_petap.pl in training); the decision is 'not enough' when
          that marker exists.

The multiplicity of a hint is the larger of column 6 and the mult=
attribute (bam2hints.gff has it in column 6, ProtHint hints in both).
"""

import argparse
import os
import re
import sys

MULT_RE = re.compile(r"(?:^|;)\s*mult=(\d+)")
INT_RE = re.compile(r"(\d+)")


def hint_multiplicity(fields):
    """Multiplicity of one GFF hint line (list of 9 columns), at least 1."""
    mult = 1
    m = INT_RE.match(fields[5].strip())
    if m:
        mult = max(mult, int(m.group(1)))
    m = MULT_RE.search(fields[8])
    if m:
        mult = max(mult, int(m.group(1)))
    return mult


def count_introns(path, min_mult):
    """Return (number of intron hints, number with multiplicity >= min_mult)."""
    n_introns = 0
    n_supported = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() != "intron":
                continue
            n_introns += 1
            if hint_multiplicity(fields) >= min_mult:
                n_supported += 1
    return n_introns, n_supported


def write_decision(path, entries):
    """Write the decision file: one key<TAB>value line per entry."""
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w") as handle:
        for key, value in entries:
            handle.write(f"{key}\t{str(value).replace(chr(9), ' ').replace(chr(10), ' ')}\n")


def read_decision(path):
    """Read a decision file into a dict (used by the tests and the report)."""
    entries = {}
    with open(path) as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition("\t")
            if key:
                entries[key] = value
    return entries


def cmd_introns(args):
    n_introns, n_supported = count_introns(args.hints, args.min_mult)
    sufficient = n_supported >= args.min_introns
    reason = (
        f"{n_supported} {args.evidence} intron hints with multiplicity >= "
        f"{args.min_mult} ({n_introns} intron hints in total); at least "
        f"{args.min_introns} are needed"
    )
    write_decision(args.out, [
        ("check", "introns"),
        ("evidence", args.evidence),
        ("sufficient", int(sufficient)),
        ("introns", n_introns),
        ("introns_supported", n_supported),
        ("min_mult", args.min_mult),
        ("min_introns", args.min_introns),
        ("reason", reason),
    ])
    print(("enough: " if sufficient else "NOT enough: ") + reason)
    return 0


def cmd_result(args):
    tool = args.tool
    if os.path.exists(args.marker):
        with open(args.marker) as handle:
            reason = " ".join(l.strip() for l in handle if l.strip())
        reason = reason or f"{tool} stopped for lack of training data"
        sufficient = False
    elif os.path.exists(args.gtf) and os.path.getsize(args.gtf) > 0:
        reason = f"{tool} finished and wrote genemark.gtf"
        sufficient = True
    else:
        sys.stderr.write(
            f"ERROR: {args.gtf} is empty or missing and there is no marker "
            f"file {args.marker}: {tool} did not finish.\n"
        )
        return 1
    write_decision(args.out, [
        ("check", "result"),
        ("evidence", tool),
        ("sufficient", int(sufficient)),
        ("reason", reason),
    ])
    print(("enough: " if sufficient else "NOT enough: ") + reason)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("introns", help="count intron hints with enough multiplicity")
    p.add_argument("--hints", required=True, help="GFF file with intron hints")
    p.add_argument("--evidence", required=True, help="label, e.g. RNA-Seq or protein")
    p.add_argument("--min-mult", type=int, required=True)
    p.add_argument("--min-introns", type=int, required=True)
    p.add_argument("--out", required=True, help="decision file")
    p.set_defaults(func=cmd_introns)

    p = sub.add_parser("result", aliases=["etp"],
                       help="decision from the result of a GeneMark run")
    p.add_argument("--gtf", required=True, help="genemark.gtf of the run")
    p.add_argument("--marker", required=True,
                   help="marker file written when the program lacked training data")
    p.add_argument("--tool", default="GeneMark-ETP",
                   help="program name for the messages (default: GeneMark-ETP)")
    p.add_argument("--out", required=True, help="decision file")
    p.set_defaults(func=cmd_result)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
