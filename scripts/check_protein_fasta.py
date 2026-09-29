#!/usr/bin/env python3
"""
Check a protein FASTA file against the rules of GeneMark-ETP's
re_format_aa_fasta.pl before GeneMark-ETP runs.

GeneMark-ETP aborts with the uninformative "error in protein file parsing:
proteins.fa" if the protein database has a duplicated ID (first word of the
header) or a sequence line with characters other than letters, '*' and '-'
(issue #99). This script reports the offending records instead and exits 1.

Usage:
    check_protein_fasta.py proteins.fa
"""

import re
import sys

MAX_EXAMPLES = 5
VALID_SEQ = re.compile(r"^[A-Z*\-]+$")


def check(path):
    seen = set()
    duplicates = []
    bad_lines = []
    orphan_lines = 0
    n_records = 0
    current = None
    with open(path, errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            if line.startswith(">"):
                m = re.match(r">\s*(\S+)", line)
                current = m.group(1) if m else ""
                n_records += 1
                if current in seen:
                    duplicates.append((lineno, current))
                seen.add(current)
                continue
            if current is None:
                orphan_lines += 1
                continue
            seq = re.sub(r"\s", "", line.upper())
            if not VALID_SEQ.match(seq):
                bad = sorted(set(re.sub(r"[A-Z*\-]", "", seq)))
                bad_lines.append((lineno, current, "".join(bad)))
    return n_records, duplicates, bad_lines, orphan_lines


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: check_protein_fasta.py proteins.fa")
    path = sys.argv[1]
    n_records, duplicates, bad_lines, orphan_lines = check(path)

    errors = []
    if n_records == 0:
        errors.append("no FASTA records found")
    if orphan_lines:
        errors.append(f"{orphan_lines} sequence line(s) before the first '>' header")
    if duplicates:
        errors.append(f"{len(duplicates)} duplicated protein ID(s) (first word of the header):")
        for lineno, pid in duplicates[:MAX_EXAMPLES]:
            errors.append(f"    line {lineno}: {pid}")
        errors.append("  Did you list the same protein file twice, or concatenate overlapping databases?")
    if bad_lines:
        errors.append(f"{len(bad_lines)} sequence line(s) with characters other than letters, '*' or '-':")
        for lineno, pid, chars in bad_lines[:MAX_EXAMPLES]:
            errors.append(f"    line {lineno} (record {pid}): {chars!r}")

    if errors:
        print(f"ERROR: protein file {path} would be rejected by GeneMark-ETP:", file=sys.stderr)
        for e in errors:
            print("  " + e, file=sys.stderr)
        sys.exit(1)
    print(f"Protein file OK: {n_records} records")


if __name__ == "__main__":
    main()
