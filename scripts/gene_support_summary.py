#!/usr/bin/env python3

"""
Compute per-gene evidence support summary from hints file.

For each gene/transcript in the GTF, counts how many introns and exons
are supported by hints from the hintsfile, split by evidence source
(RNA-Seq vs protein).

Hint source mapping:
  E = RNA-Seq (bam2hints, intron hints)
  P = protein (ProtHint)
  M = manual/compleasm
  W = wiggle (coverage)
  C = combined (ETP)

Output: TSV with per-transcript support statistics.

Usage:
    python3 gene_support_summary.py -g braker.gtf -H hintsfile.gff -o gene_support.tsv
"""

import argparse
import sys
import re
from collections import defaultdict
from intervaltree import IntervalTree


def parse_hints(hints_file):
    """Parse hints file into interval trees by (chrom, feature_type, source_class).

    Returns dict: (chrom,) -> IntervalTree of (start, end, source_class)
    where source_class is 'rnaseq', 'protein', or 'other'.
    """
    # Intron hints indexed by chrom
    intron_hints = defaultdict(IntervalTree)
    # Exon-like hints (exonpart, CDSpart, ep) indexed by chrom
    exon_hints = defaultdict(IntervalTree)

    src_pattern = re.compile(r'src=([^;]+)')

    rnaseq_sources = {'E', 'W', 'b2h'}
    protein_sources = {'P', 'PH'}
    combined_sources = {'C'}  # ETP combined hints — count as both

    with open(hints_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]

            # Determine source class
            src_match = src_pattern.search(attrs)
            src = src_match.group(1) if src_match else ''

            if src in combined_sources:
                source_class = 'combined'
            elif src in rnaseq_sources:
                source_class = 'rnaseq'
            elif src in protein_sources:
                source_class = 'protein'
            else:
                source_class = 'other'

            if feature == 'intron':
                # Store as exact interval for intron matching
                intron_hints[chrom].addi(start, end + 1, source_class)
            elif feature in ('exonpart', 'CDSpart', 'ep', 'exon', 'CDS'):
                exon_hints[chrom].addi(start, end + 1, source_class)

    return intron_hints, exon_hints


def parse_gtf_transcripts(gtf_file):
    """Parse GTF into per-transcript feature lists.

    Returns list of (gene_id, tx_id, chrom, strand, introns, exons)
    where introns and exons are lists of (start, end).
    """
    tx_id_pattern = re.compile(r'transcript_id "([^"]+)"')
    gene_id_pattern = re.compile(r'gene_id "([^"]+)"')

    # Collect CDS/exon features per transcript
    tx_features = defaultdict(list)  # tx_id -> [(feature, start, end, chrom, strand)]
    tx_gene = {}  # tx_id -> gene_id

    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature not in ('CDS', 'exon'):
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]

            tx_match = tx_id_pattern.search(attrs)
            gene_match = gene_id_pattern.search(attrs)
            if not tx_match:
                continue

            tx_id = tx_match.group(1)
            gene_id = gene_match.group(1) if gene_match else 'unknown'
            tx_gene[tx_id] = gene_id
            tx_features[tx_id].append((feature, start, end, chrom, strand))

    # Derive introns from exon coordinates
    results = []
    for tx_id, features in tx_features.items():
        exons = sorted([(s, e) for f, s, e, c, st in features if f == 'exon'],
                       key=lambda x: x[0])
        cds = sorted([(s, e) for f, s, e, c, st in features if f == 'CDS'],
                     key=lambda x: x[0])

        if not exons and not cds:
            continue

        chrom = features[0][3]
        strand = features[0][4]

        # Derive introns from exons
        introns = []
        ref = exons if exons else cds
        for i in range(len(ref) - 1):
            intron_start = ref[i][1] + 1
            intron_end = ref[i + 1][0] - 1
            if intron_start <= intron_end:
                introns.append((intron_start, intron_end))

        results.append((tx_gene[tx_id], tx_id, chrom, strand, introns, exons))

    return results


def check_support(intervals, hints_tree, chrom):
    """Check how many intervals are supported by hints.

    For introns: exact match (start and end must match a hint).
    Returns (total, supported_rnaseq, supported_protein, supported_any).
    """
    total = len(intervals)
    sup_rnaseq = 0
    sup_protein = 0
    sup_any = 0

    tree = hints_tree.get(chrom, IntervalTree())

    for start, end in intervals:
        # Query for overlapping hints at this exact position
        overlaps = tree.overlap(start, end + 1)

        has_rnaseq = False
        has_protein = False

        for iv in overlaps:
            # For introns, require exact match
            if iv.begin == start and iv.end == end + 1:
                if iv.data == 'rnaseq':
                    has_rnaseq = True
                elif iv.data == 'protein':
                    has_protein = True
                elif iv.data == 'combined':
                    has_rnaseq = True
                    has_protein = True

        if has_rnaseq:
            sup_rnaseq += 1
        if has_protein:
            sup_protein += 1
        if has_rnaseq or has_protein:
            sup_any += 1

    return total, sup_rnaseq, sup_protein, sup_any


def check_exon_support(exons, hints_tree, chrom):
    """Check how many exons have any overlapping exon/CDS hints.

    Uses overlap (not exact match) since exonpart hints cover partial exons.
    Returns (total, supported_rnaseq, supported_protein, supported_any).
    """
    total = len(exons)
    sup_rnaseq = 0
    sup_protein = 0
    sup_any = 0

    tree = hints_tree.get(chrom, IntervalTree())

    for start, end in exons:
        overlaps = tree.overlap(start, end + 1)

        has_rnaseq = any(iv.data in ('rnaseq', 'combined') for iv in overlaps)
        has_protein = any(iv.data in ('protein', 'combined') for iv in overlaps)

        if has_rnaseq:
            sup_rnaseq += 1
        if has_protein:
            sup_protein += 1
        if has_rnaseq or has_protein:
            sup_any += 1

    return total, sup_rnaseq, sup_protein, sup_any


def main():
    parser = argparse.ArgumentParser(
        description="Compute per-gene evidence support summary.")
    parser.add_argument("-g", "--gtf", required=True,
                        help="Input GTF file (braker.gtf)")
    parser.add_argument("-H", "--hints", required=True,
                        help="Hints file (hintsfile.gff)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file")
    args = parser.parse_args()

    print("Loading hints...", file=sys.stderr)
    intron_hints, exon_hints = parse_hints(args.hints)

    print("Parsing gene structures...", file=sys.stderr)
    transcripts = parse_gtf_transcripts(args.gtf)

    print(f"Computing support for {len(transcripts)} transcripts...", file=sys.stderr)

    # Collect per-transcript rows and accumulate totals
    rows = []
    tot_introns = 0
    tot_i_rna = 0
    tot_i_prot = 0
    tot_i_any = 0
    tot_exons = 0
    tot_e_rna = 0
    tot_e_prot = 0
    tot_e_any = 0
    # Transcript-level counts
    n_tx = 0
    n_tx_fully_supported = 0       # all introns supported by any evidence
    n_tx_any_intron_supported = 0  # at least one intron supported
    n_tx_no_intron_support = 0     # multi-exon but zero intron support
    n_tx_single_exon = 0           # no introns to evaluate

    for gene_id, tx_id, chrom, strand, introns, exons in transcripts:
        n_introns, sup_i_rna, sup_i_prot, sup_i_any = check_support(
            introns, intron_hints, chrom)
        n_exons, sup_e_rna, sup_e_prot, sup_e_any = check_exon_support(
            exons, exon_hints, chrom)

        frac = f"{sup_i_any / n_introns:.3f}" if n_introns > 0 else "NA"

        pct_i_rna = f"{100 * sup_i_rna / n_introns:.1f}" if n_introns > 0 else "NA"
        pct_i_prot = f"{100 * sup_i_prot / n_introns:.1f}" if n_introns > 0 else "NA"
        pct_i_any = f"{100 * sup_i_any / n_introns:.1f}" if n_introns > 0 else "NA"
        pct_e_rna = f"{100 * sup_e_rna / n_exons:.1f}" if n_exons > 0 else "NA"
        pct_e_prot = f"{100 * sup_e_prot / n_exons:.1f}" if n_exons > 0 else "NA"
        pct_e_any = f"{100 * sup_e_any / n_exons:.1f}" if n_exons > 0 else "NA"

        rows.append((gene_id, tx_id, chrom, strand,
                      n_introns, sup_i_rna, pct_i_rna, sup_i_prot, pct_i_prot, sup_i_any, pct_i_any,
                      n_exons, sup_e_rna, pct_e_rna, sup_e_prot, pct_e_prot, sup_e_any, pct_e_any))

        tot_introns += n_introns
        tot_i_rna += sup_i_rna
        tot_i_prot += sup_i_prot
        tot_i_any += sup_i_any
        tot_exons += n_exons
        tot_e_rna += sup_e_rna
        tot_e_prot += sup_e_prot
        tot_e_any += sup_e_any
        n_tx += 1
        if n_introns == 0:
            n_tx_single_exon += 1
        elif sup_i_any == n_introns:
            n_tx_fully_supported += 1
        elif sup_i_any > 0:
            n_tx_any_intron_supported += 1
        else:
            n_tx_no_intron_support += 1

    def pct(num, denom):
        return f"{100 * num / denom:.1f}%" if denom > 0 else "NA"

    with open(args.output, 'w') as out:
        # Summary header
        out.write(f"# Gene support summary\n")
        out.write(f"# Transcripts: {n_tx}\n")
        out.write(f"#   Single-exon (no introns to evaluate): {n_tx_single_exon} ({pct(n_tx_single_exon, n_tx)})\n")
        out.write(f"#   All introns supported: {n_tx_fully_supported} ({pct(n_tx_fully_supported, n_tx)})\n")
        out.write(f"#   Some introns supported: {n_tx_any_intron_supported} ({pct(n_tx_any_intron_supported, n_tx)})\n")
        out.write(f"#   No intron support: {n_tx_no_intron_support} ({pct(n_tx_no_intron_support, n_tx)})\n")
        out.write(f"# Introns: {tot_introns}\n")
        out.write(f"#   Supported by RNA-Seq: {tot_i_rna} ({pct(tot_i_rna, tot_introns)})\n")
        out.write(f"#   Supported by protein: {tot_i_prot} ({pct(tot_i_prot, tot_introns)})\n")
        out.write(f"#   Supported by any: {tot_i_any} ({pct(tot_i_any, tot_introns)})\n")
        out.write(f"# Exons: {tot_exons}\n")
        out.write(f"#   Supported by RNA-Seq: {tot_e_rna} ({pct(tot_e_rna, tot_exons)})\n")
        out.write(f"#   Supported by protein: {tot_e_prot} ({pct(tot_e_prot, tot_exons)})\n")
        out.write(f"#   Supported by any: {tot_e_any} ({pct(tot_e_any, tot_exons)})\n")

        # Column header
        out.write("gene_id\ttranscript_id\tchrom\tstrand\t"
                  "num_introns\t"
                  "introns_sup_rnaseq\tpct_introns_sup_rnaseq\t"
                  "introns_sup_protein\tpct_introns_sup_protein\t"
                  "introns_sup_any\tpct_introns_sup_any\t"
                  "num_exons\t"
                  "exons_sup_rnaseq\tpct_exons_sup_rnaseq\t"
                  "exons_sup_protein\tpct_exons_sup_protein\t"
                  "exons_sup_any\tpct_exons_sup_any\n")

        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
