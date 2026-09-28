#!/usr/bin/env python3
"""
Evaluate a gene set against the held-out AUGUSTUS test set (train.gb.test).

The test set is the portion of the training genes that split_training_set
keeps away from etraining and optimize_augustus.pl. AUGUSTUS itself is scored
on it by predicting ab initio on each GenBank sequence. This script scores an
arbitrary GTF (normally the final braker.gtf) on the same loci, so the numbers
sit next to accuracy_after_training.txt / accuracy_after_optimize.txt.

How it works:
  - Each GenBank entry is one genomic region. gff2gbSmallDNA.pl names it
    "<seqname>_<start>-<end>" (1-based, inclusive), so the region and the
    annotated CDS structures can be lifted back to genome coordinates.
  - Predicted transcripts are restricted to those regions. CDS parts outside
    a region are clipped, as they would be if AUGUSTUS predicted on the
    GenBank sequence alone.
  - Nucleotide level: coding bases (strand aware).
  - Exon level: CDS exons with identical start, end and strand.
  - Gene level: an annotated gene is found if some predicted transcript has
    exactly its CDS exon chain; a predicted gene is correct if any of its
    transcripts matches an annotated gene.

This is a QC number, not an independent benchmark: the test genes come from
GeneMark (or GeneMark-ETP training genes), and the final gene set contains
GeneMark transcripts selected by TSEBRA, so agreement is biased upwards.

The output follows the format of accuracy_after_optimize.txt so that
training_summary.py can parse it.

Usage:
    eval_test_set_accuracy.py -t train.gb.test -g braker.gtf -o accuracy_final_gene_set.txt
"""

import argparse
import re
import sys
from collections import defaultdict

LOCUS_RE = re.compile(r'^LOCUS\s+(\S+)_(\d+)-(\d+)\s+(\d+)\s+bp')


def parse_location(loc):
    """Parse a GenBank location string into (strand, [(start, end), ...]).

    Handles join(...), complement(...) and complement(join(...)). Partial
    markers < and > are ignored. Coordinates are 1-based and relative to the
    GenBank sequence.
    """
    loc = loc.replace(' ', '').replace('<', '').replace('>', '')
    strand = '+'
    if loc.startswith('complement(') and loc.endswith(')'):
        strand = '-'
        loc = loc[len('complement('):-1]
    if loc.startswith('join(') and loc.endswith(')'):
        loc = loc[len('join('):-1]
    parts = []
    for piece in loc.split(','):
        if '..' in piece:
            a, b = piece.split('..')
        else:
            a = b = piece
        parts.append((int(a), int(b)))
    return strand, sorted(parts)


def parse_genbank_test_set(path):
    """Return a list of regions from a gff2gbSmallDNA.pl GenBank file.

    Each region is a dict with seq, start, end and genes, where genes is a
    list of (strand, tuple_of_cds_exons) in genome coordinates.
    """
    regions = []
    region = None
    in_features = False
    feature = None      # (key, location_string) of the feature being read
    loc_done = False    # qualifiers have started; the location is complete

    def finish_feature():
        if region is None or feature is None:
            return
        key, loc = feature
        if key != 'CDS':
            return
        strand, parts = parse_location(loc)
        offset = region['start'] - 1
        exons = tuple((s + offset, e + offset) for s, e in parts)
        region['genes'].append((strand, exons))

    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('LOCUS'):
                m = LOCUS_RE.match(line)
                if not m:
                    raise ValueError(f"Cannot parse genome coordinates from LOCUS line: {line}")
                seq, start, end, length = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
                if end - start + 1 != length:
                    raise ValueError(
                        f"LOCUS {seq}_{start}-{end} has length {length}, expected {end - start + 1}"
                    )
                region = {'seq': seq, 'start': start, 'end': end, 'genes': []}
                regions.append(region)
                in_features = False
                feature = None
            elif line.startswith('FEATURES'):
                in_features = True
            elif line.startswith('BASE COUNT') or line.startswith('ORIGIN') or line.startswith('//'):
                finish_feature()
                feature = None
                in_features = False
            elif in_features:
                key = line[5:21].strip()
                value = line[21:].strip()
                if key:
                    finish_feature()
                    feature = (key, value)
                    loc_done = False
                elif value.startswith('/'):
                    loc_done = True
                elif feature is not None and not loc_done:
                    # Location continued on the next line
                    feature = (feature[0], feature[1] + value)
    return regions


def parse_gtf_cds(path):
    """Return {seq: [(gene_id, transcript_id, strand, [(start, end), ...])]}.

    GenBank CDS features from gff2gbSmallDNA.pl include the stop codon, while
    AUGUSTUS-style GTF keeps it as a separate stop_codon feature. Stop codons
    are therefore merged into the CDS chain.
    """
    tx = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] not in ('CDS', 'stop_codon'):
                continue
            m_tx = re.search(r'transcript_id "([^"]+)"', cols[8])
            m_g = re.search(r'gene_id "([^"]+)"', cols[8])
            if not m_tx:
                continue
            tid = m_tx.group(1)
            gid = m_g.group(1) if m_g else tid
            key = (cols[0], gid, tid, cols[6])
            tx.setdefault(key, []).append((int(cols[3]), int(cols[4])))
    by_seq = defaultdict(list)
    for (seq, gid, tid, strand), cds in tx.items():
        exons = [tuple(iv) for iv in merge_intervals(cds)]
        by_seq[seq].append((gid, tid, strand, exons))
    return by_seq


def merge_intervals(intervals):
    """Merge overlapping or adjacent intervals; return total length and merged list."""
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def overlap_length(a, b):
    """Total overlap of two merged, sorted interval lists."""
    i = j = total = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s <= e:
            total += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def length(intervals):
    return sum(e - s + 1 for s, e in intervals)


def evaluate(regions, predictions):
    """Compute counts for nucleotide, exon and gene level accuracy."""
    c = defaultdict(int)
    c['regions'] = len(regions)
    for r in regions:
        rs, re_ = r['start'], r['end']
        # Clip predicted transcripts to the region
        pred_tx = []
        for gid, tid, strand, cds in predictions.get(r['seq'], []):
            clipped = tuple((max(s, rs), min(e, re_)) for s, e in cds if e >= rs and s <= re_)
            if clipped:
                pred_tx.append((gid, strand, clipped))

        for strand in '+-':
            anno_bases = merge_intervals([x for st, ex in r['genes'] if st == strand for x in ex])
            pred_bases = merge_intervals([x for _, st, ex in pred_tx if st == strand for x in ex])
            c['nu_tp'] += overlap_length(anno_bases, pred_bases)
            c['nu_anno'] += length(anno_bases)
            c['nu_pred'] += length(pred_bases)

        anno_exons = {(st, s, e) for st, ex in r['genes'] for s, e in ex}
        pred_exons = {(st, s, e) for _, st, ex in pred_tx for s, e in ex}
        c['ex_tp'] += len(anno_exons & pred_exons)
        c['ex_anno'] += len(anno_exons)
        c['ex_pred'] += len(pred_exons)

        anno_chains = {(st, ex) for st, ex in r['genes']}
        pred_chains = {(st, ex) for _, st, ex in pred_tx}
        c['gene_anno'] += len(anno_chains)
        c['gene_tp_anno'] += len(anno_chains & pred_chains)
        pred_genes = defaultdict(set)
        for gid, st, ex in pred_tx:
            pred_genes[gid].add((st, ex))
        c['gene_pred'] += len(pred_genes)
        c['gene_tp_pred'] += sum(1 for chains in pred_genes.values() if chains & anno_chains)
    return c


def pct(num, den):
    return 100.0 * num / den if den else 0.0


def format_report(c, gtf_name):
    nu_sen, nu_sp = pct(c['nu_tp'], c['nu_anno']), pct(c['nu_tp'], c['nu_pred'])
    ex_sen, ex_sp = pct(c['ex_tp'], c['ex_anno']), pct(c['ex_tp'], c['ex_pred'])
    gen_sen, gen_sp = pct(c['gene_tp_anno'], c['gene_anno']), pct(c['gene_tp_pred'], c['gene_pred'])
    target = (3 * nu_sen + 2 * nu_sp + 4 * ex_sen + 3 * ex_sp + 2 * gen_sen + 1 * gen_sp) / 15
    return f"""Accuracy of the final gene set on the held-out AUGUSTUS test set

Gene set: {gtf_name}
Test loci: {c['regions']}
Test genes: {c['gene_anno']}
Predicted genes in test loci: {c['gene_pred']}
Test genes reproduced exactly (CDS): {c['gene_tp_anno']}

Nucleotide level:
  Sensitivity: {nu_sen:.2f}%
  Specificity: {nu_sp:.2f}%

Exon level:
  Sensitivity: {ex_sen:.2f}%
  Specificity: {ex_sp:.2f}%

Gene level:
  Sensitivity: {gen_sen:.2f}%
  Specificity: {gen_sp:.2f}%

Target accuracy (weighted): {target:.2f}%
Formula: (3*nu_sen + 2*nu_sp + 4*ex_sen + 3*ex_sp + 2*gen_sen + 1*gen_sp) / 15

Note: QC only, not an independent benchmark. The test genes are GeneMark
training genes and the final gene set contains GeneMark transcripts chosen
by TSEBRA, so these values are biased upwards compared to the AUGUSTUS-only
values in accuracy_after_optimize.txt.
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('-t', '--test-gb', required=True, help='Held-out GenBank test set (train.gb.test)')
    ap.add_argument('-g', '--gtf', required=True, help='Gene set to evaluate (GTF with CDS features)')
    ap.add_argument('-o', '--output', required=True, help='Output accuracy file')
    args = ap.parse_args()

    regions = parse_genbank_test_set(args.test_gb)
    if not regions:
        sys.exit(f"[ERROR] No LOCUS entries in {args.test_gb}")
    predictions = parse_gtf_cds(args.gtf)
    counts = evaluate(regions, predictions)
    report = format_report(counts, args.gtf)
    with open(args.output, 'w') as f:
        f.write(report)
    print(report)


if __name__ == '__main__':
    main()
