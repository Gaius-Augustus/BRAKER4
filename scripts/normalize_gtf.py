#!/usr/bin/env python3

"""
Normalize CDS boundaries and validate gene structures in a GTF file.

Fixes:
- Stop codon included in CDS (GeneMark convention) — trims 3bp from terminal CDS
- Validates all gene structures after trimming; discards broken genes
- Reports transcripts with non-ATG start codons (BRAKER issue #283)
- Reports transcripts with CDS not divisible by 3

Usage:
    python3 normalize_gtf.py -g genome.fa -f braker.gtf -o braker.normalized.gtf -l normalize.log

Author: Generated for BRAKER4 pipeline
"""

import argparse
import sys
from collections import defaultdict
from Bio import SeqIO


def parse_gtf(gtf_file):
    """Parse GTF into gene -> transcript -> features structure."""
    genes = defaultdict(lambda: defaultdict(list))  # gene_id -> tx_id -> [features]
    gene_lines = {}  # gene_id -> gene line
    tx_lines = {}    # tx_id -> transcript line
    tx_to_gene = {}  # tx_id -> gene_id

    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get('gene_id')
            tx_id = attrs.get('transcript_id')

            if feature == 'gene':
                gene_lines[gene_id or fields[8].strip()] = fields
            elif feature == 'transcript':
                key = tx_id or fields[8].strip()
                tx_lines[key] = fields
                if gene_id:
                    tx_to_gene[key] = gene_id
            elif tx_id and gene_id:
                genes[gene_id][tx_id].append(fields)
                tx_to_gene[tx_id] = gene_id

    return genes, gene_lines, tx_lines, tx_to_gene


def parse_attrs(attr_str):
    """Parse GTF attribute string into dict."""
    attrs = {}
    for attr in attr_str.split(';'):
        attr = attr.strip()
        if not attr:
            continue
        parts = attr.split(' ', 1)
        if len(parts) == 2:
            attrs[parts[0]] = parts[1].strip('"')
    return attrs


def get_sequence(genome, chrom, start, end, strand):
    """Extract sequence from genome (1-based coords)."""
    if chrom not in genome:
        return None
    seq = genome[chrom].seq[start - 1:end]
    if strand == '-':
        seq = seq.reverse_complement()
    return str(seq).upper()


def get_cds_features(features):
    """Get CDS features sorted by genomic position."""
    cds = [f for f in features if f[2] == 'CDS']
    cds.sort(key=lambda f: int(f[3]))
    return cds


def normalize_transcript(features, genome, log_messages):
    """
    Normalize CDS boundaries for a transcript.

    Returns (normalized_features, is_valid, messages).
    """
    cds_features = get_cds_features(features)
    if not cds_features:
        return features, True, []

    strand = cds_features[0][6]
    tx_id = parse_attrs(cds_features[0][8]).get('transcript_id', '?')
    chrom = cds_features[0][0]
    messages = []

    # Determine terminal CDS (last in translation order)
    if strand == '+':
        terminal_cds = cds_features[-1]
    else:
        terminal_cds = cds_features[0]

    terminal_start = int(terminal_cds[3])
    terminal_end = int(terminal_cds[4])

    # Check if terminal CDS ends in a stop codon
    if strand == '+':
        stop_seq = get_sequence(genome, chrom, terminal_end - 2, terminal_end, '+')
    else:
        stop_seq = get_sequence(genome, chrom, terminal_start, terminal_start + 2, '-')

    stop_codons = {'TAA', 'TAG', 'TGA'}

    if stop_seq and stop_seq in stop_codons:
        # Trim stop codon from CDS
        if strand == '+':
            new_end = terminal_end - 3
            if new_end < terminal_start:
                messages.append(f"DISCARD {tx_id}: trimming stop codon would make CDS length <= 0")
                return features, False, messages
            terminal_cds[4] = str(new_end)
        else:
            new_start = terminal_start + 3
            if new_start > terminal_end:
                messages.append(f"DISCARD {tx_id}: trimming stop codon would make CDS length <= 0")
                return features, False, messages
            terminal_cds[3] = str(new_start)
        messages.append(f"TRIMMED {tx_id}: removed stop codon ({stop_seq}) from CDS")

    # Recalculate CDS features after trimming
    cds_features = get_cds_features([f for f in features if f[2] == 'CDS'])

    # Validate: all CDS must have start <= end
    for cds in cds_features:
        if int(cds[3]) > int(cds[4]):
            messages.append(f"DISCARD {tx_id}: CDS start > end after trimming ({cds[3]} > {cds[4]})")
            return features, False, messages

    # Validate: CDS must be contained within exons
    exons = [f for f in features if f[2] == 'exon']
    for cds in cds_features:
        cds_s, cds_e = int(cds[3]), int(cds[4])
        contained = any(int(e[3]) <= cds_s and cds_e <= int(e[4]) for e in exons)
        if not contained:
            messages.append(f"DISCARD {tx_id}: CDS ({cds_s}-{cds_e}) not contained in any exon")
            return features, False, messages

    # Validate: total CDS length divisible by 3
    total_cds_len = sum(int(c[4]) - int(c[3]) + 1 for c in cds_features)
    if total_cds_len % 3 != 0:
        messages.append(f"WARNING {tx_id}: total CDS length {total_cds_len} not divisible by 3")
        # Don't discard — this can happen with partial genes

    # Check start codon (issue #283)
    if strand == '+':
        first_cds = cds_features[0]
        start_seq = get_sequence(genome, chrom, int(first_cds[3]), int(first_cds[3]) + 2, '+')
    else:
        first_cds = cds_features[-1]
        start_seq = get_sequence(genome, chrom, int(first_cds[4]) - 2, int(first_cds[4]), '-')

    if start_seq and start_seq != 'ATG':
        messages.append(f"WARNING {tx_id}: non-ATG start codon ({start_seq})")

    return features, True, messages


def update_gene_boundaries(gene_id, genes, gene_lines, tx_lines):
    """Update gene and transcript boundaries to match their features."""
    all_starts = []
    all_ends = []

    for tx_id, features in genes[gene_id].items():
        tx_starts = [int(f[3]) for f in features]
        tx_ends = [int(f[4]) for f in features]
        if tx_starts and tx_ends:
            tx_min = min(tx_starts)
            tx_max = max(tx_ends)
            all_starts.append(tx_min)
            all_ends.append(tx_max)
            if tx_id in tx_lines:
                tx_lines[tx_id][3] = str(tx_min)
                tx_lines[tx_id][4] = str(tx_max)

    if all_starts and all_ends and gene_id in gene_lines:
        gene_lines[gene_id][3] = str(min(all_starts))
        gene_lines[gene_id][4] = str(max(all_ends))


def write_gtf(output_file, genes, gene_lines, tx_lines, tx_to_gene):
    """Write normalized GTF, preserving gene/transcript/feature order."""
    # Group transcripts by gene
    gene_order = []
    seen_genes = set()

    # Determine gene order by genomic position
    for gene_id in sorted(gene_lines.keys(),
                          key=lambda g: (gene_lines[g][0], int(gene_lines[g][3]))):
        if gene_id in genes:
            gene_order.append(gene_id)

    with open(output_file, 'w') as f:
        for gene_id in gene_order:
            if gene_id not in genes or not genes[gene_id]:
                continue

            # Write gene line
            if gene_id in gene_lines:
                f.write('\t'.join(gene_lines[gene_id]) + '\n')

            # Write transcripts
            for tx_id in sorted(genes[gene_id].keys()):
                if tx_id in tx_lines:
                    f.write('\t'.join(tx_lines[tx_id]) + '\n')
                for feature in sorted(genes[gene_id][tx_id],
                                      key=lambda x: (int(x[3]), x[2])):
                    f.write('\t'.join(feature) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Normalize CDS boundaries and validate gene structures.")
    parser.add_argument("-g", "--genome", required=True,
                        help="Genome FASTA file")
    parser.add_argument("-f", "--gtf", required=True,
                        help="Input GTF file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output normalized GTF file")
    parser.add_argument("-l", "--log", required=True,
                        help="Log file for normalization report")
    args = parser.parse_args()

    # Load genome
    print(f"Loading genome: {args.genome}", file=sys.stderr)
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    # Parse GTF
    print(f"Parsing GTF: {args.gtf}", file=sys.stderr)
    genes, gene_lines, tx_lines, tx_to_gene = parse_gtf(args.gtf)

    # Normalize each transcript
    all_messages = []
    discarded_genes = set()
    trimmed_count = 0
    non_atg_count = 0
    frame_warn_count = 0

    for gene_id, transcripts in list(genes.items()):
        gene_valid = True
        for tx_id, features in list(transcripts.items()):
            features, is_valid, messages = normalize_transcript(
                features, genome, all_messages)
            all_messages.extend(messages)

            for msg in messages:
                if msg.startswith("TRIMMED"):
                    trimmed_count += 1
                elif msg.startswith("WARNING") and "non-ATG" in msg:
                    non_atg_count += 1
                elif msg.startswith("WARNING") and "divisible" in msg:
                    frame_warn_count += 1

            if not is_valid:
                gene_valid = False
                break

        if not gene_valid:
            discarded_genes.add(gene_id)
            del genes[gene_id]
        else:
            update_gene_boundaries(gene_id, genes, gene_lines, tx_lines)

    # Write output
    write_gtf(args.output, genes, gene_lines, tx_lines, tx_to_gene)

    # Write log
    genes_after = len(genes)
    with open(args.log, 'w') as log:
        log.write(f"CDS boundary normalization report\n")
        log.write(f"{'=' * 50}\n")
        log.write(f"Stop codons trimmed from CDS: {trimmed_count}\n")
        log.write(f"Genes discarded (broken after trim): {len(discarded_genes)}\n")
        log.write(f"Non-ATG start codons (warning only): {non_atg_count}\n")
        log.write(f"CDS length not divisible by 3 (warning): {frame_warn_count}\n")
        log.write(f"Genes in output: {genes_after}\n")
        log.write(f"\nDetails:\n")
        for msg in all_messages:
            log.write(f"  {msg}\n")

    print(f"Done: {trimmed_count} trimmed, {len(discarded_genes)} discarded, "
          f"{non_atg_count} non-ATG starts, {genes_after} genes in output",
          file=sys.stderr)


if __name__ == "__main__":
    main()
