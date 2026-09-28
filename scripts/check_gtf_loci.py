#!/usr/bin/env python3

"""
Check that every gene model in a GTF sits on a single locus (issue #62).

Fails (exit code 1) when
- a transcript has features on more than one sequence or strand,
- a transcript line differs from its features in sequence, strand or span,
- a gene line differs from its transcripts in sequence, strand or span,
- a gene or transcript line appears more than once.

Runs on braker.gtf and braker_utr.gtf before GFF3 conversion, so a broken
model stops the pipeline instead of reaching braker.gff3.

Handles attribute-style column 9 as well as the bare TSEBRA/AUGUSTUS gene
and transcript lines, where column 9 holds only the ID.

Usage:
    python3 check_gtf_loci.py braker.gtf
"""

import argparse
import re
import sys
from collections import defaultdict

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')
MAX_REPORTED = 20


def _id(attr, pattern):
    m = pattern.search(attr)
    return m.group(1) if m else attr.strip()


def _attr(attr, pattern):
    m = pattern.search(attr)
    return m.group(1) if m else None


def parse(gtf_file):
    """Return gene lines, transcript lines, features per transcript and
    transcript -> gene links. Lines are (seq, strand, start, end)."""
    gene_lines = defaultdict(list)
    tx_lines = defaultdict(list)
    features = defaultdict(list)
    tx_gene = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            loc = (fields[0], fields[6], int(fields[3]), int(fields[4]))
            feature = fields[2]
            if feature == 'gene':
                gene_lines[_id(fields[8], GENE_RE)].append(loc)
            elif feature in ('transcript', 'mRNA'):
                tx_id = _id(fields[8], TX_RE)
                tx_lines[tx_id].append(loc)
                gene_id = _attr(fields[8], GENE_RE)
                if gene_id:
                    tx_gene[tx_id] = gene_id
            else:
                tx_id = _attr(fields[8], TX_RE)
                if tx_id is None:
                    continue
                features[tx_id].append(loc)
                gene_id = _attr(fields[8], GENE_RE)
                if gene_id:
                    tx_gene.setdefault(tx_id, gene_id)
    return gene_lines, tx_lines, features, tx_gene


def _span(locs):
    return min(l[2] for l in locs), max(l[3] for l in locs)


def _fmt(seq, strand, start, end):
    return f"{seq}:{start}-{end}({strand})"


def check(gene_lines, tx_lines, features, tx_gene):
    """Return a list of violation messages."""
    errors = []

    for kind, lines in (("gene", gene_lines), ("transcript", tx_lines)):
        for id_, locs in lines.items():
            if len(locs) > 1:
                errors.append(f"{kind} {id_}: {len(locs)} {kind} lines "
                              f"({', '.join(_fmt(*l) for l in locs)})")

    # Transcript locus: from its features, else from its transcript line
    tx_locus = {}
    for tx_id in set(features) | set(tx_lines):
        feats = features.get(tx_id, [])
        loci = sorted({(l[0], l[1]) for l in feats})
        if len(loci) > 1:
            errors.append(f"transcript {tx_id}: features on several sequences "
                          f"or strands ({', '.join(f'{s}({st})' for s, st in loci)})")
            continue
        line = tx_lines.get(tx_id, [None])[0]
        if feats:
            start, end = _span(feats)
            locus = (loci[0][0], loci[0][1], start, end)
            if line and line != locus:
                errors.append(f"transcript {tx_id}: line {_fmt(*line)} does not "
                              f"match its features {_fmt(*locus)}")
        else:
            locus = line
        tx_locus[tx_id] = locus

    gene_tx = defaultdict(list)
    for tx_id, gene_id in tx_gene.items():
        gene_tx[gene_id].append(tx_id)

    for gene_id, locs in gene_lines.items():
        txs = [tx_locus[t] for t in gene_tx.get(gene_id, []) if t in tx_locus]
        if not txs:
            continue
        loci = sorted({(l[0], l[1]) for l in txs})
        if len(loci) > 1:
            errors.append(f"gene {gene_id}: transcripts on several sequences or "
                          f"strands ({', '.join(f'{s}({st})' for s, st in loci)})")
            continue
        start, end = _span(txs)
        expected = (loci[0][0], loci[0][1], start, end)
        if locs[0] != expected:
            errors.append(f"gene {gene_id}: line {_fmt(*locs[0])} does not match "
                          f"its transcripts {_fmt(*expected)}")

    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Fail if a gene or transcript in a GTF spans more than one locus.")
    parser.add_argument("gtf", help="GTF file to check")
    args = parser.parse_args()

    gene_lines, tx_lines, features, tx_gene = parse(args.gtf)
    errors = check(gene_lines, tx_lines, features, tx_gene)
    n_tx = len(set(features) | set(tx_lines))

    if errors:
        print(f"ERROR: {args.gtf}: {len(errors)} gene models are not on a single "
              f"locus (see github.com/Gaius-Augustus/BRAKER4/issues/62):",
              file=sys.stderr)
        for msg in errors[:MAX_REPORTED]:
            print(f"  {msg}", file=sys.stderr)
        if len(errors) > MAX_REPORTED:
            print(f"  ... and {len(errors) - MAX_REPORTED} more", file=sys.stderr)
        sys.exit(1)

    print(f"Locus check passed: {len(gene_lines)} genes, {n_tx} transcripts "
          f"in {args.gtf}", file=sys.stderr)


if __name__ == "__main__":
    main()
