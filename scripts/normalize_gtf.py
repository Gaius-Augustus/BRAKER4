#!/usr/bin/env python3

"""
Normalize CDS boundaries and validate gene structures in a GTF file.

Fixes:
- Stop codon included in CDS (GeneMark convention) — trims 3bp from terminal CDS
- Validates all gene structures after trimming; discards broken genes
- Reports transcripts with non-ATG start codons (BRAKER issue #283)
- Reports transcripts with CDS not divisible by 3
- Splits gene IDs used on several loci into one gene per locus (#97)
- Removes transcripts whose CDS chain repeats an earlier transcript (#71)

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
            elif feature in ('transcript', 'mRNA'):
                key = tx_id or fields[8].strip()
                if key not in tx_lines:
                    # Normalize mRNA feature type to transcript so braker.gtf
                    # uses a single consistent feature type throughout.
                    normalized = list(fields)
                    normalized[2] = 'transcript'
                    tx_lines[key] = normalized
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


def rename_ids(fields, old_gene, new_gene, old_tx=None, new_tx=None):
    """Return a copy of a GTF line with gene and transcript IDs replaced.

    Handles attribute-style column 9 as well as the bare TSEBRA/AUGUSTUS
    gene and transcript lines, where column 9 holds only the ID.
    """
    out = list(fields)
    attr = out[8]
    if 'gene_id "' in attr or 'transcript_id "' in attr:
        attr = attr.replace(f'gene_id "{old_gene}"', f'gene_id "{new_gene}"')
        if old_tx is not None:
            attr = attr.replace(f'transcript_id "{old_tx}"',
                                f'transcript_id "{new_tx}"')
    elif old_tx is not None and attr.strip() == old_tx:
        attr = new_tx
    elif attr.strip() == old_gene:
        attr = new_gene
    out[8] = attr
    return out


def split_multilocus_ids(genes, gene_lines, tx_lines, tx_to_gene):
    """Give every locus its own gene ID (issue #97).

    A transcript ID whose features sit on more than one sequence or strand is
    split into one transcript per sequence and strand. A gene ID whose
    transcripts sit on more than one sequence or strand, or do not overlap,
    is split into one gene per locus: the first locus in genomic order keeps
    the ID, the others become <gene_id>_2, <gene_id>_3, ... and their
    transcripts are renamed to match.

    Returns log messages.
    """
    messages = []
    used_genes = set(genes) | set(gene_lines)
    used_tx = set(tx_to_gene)

    def free_id(candidate, used):
        new_id, n = candidate, 1
        while new_id in used:
            n += 1
            new_id = f"{candidate}_{n}"
        used.add(new_id)
        return new_id

    for gene_id in list(genes):
        # One piece per (transcript, sequence, strand)
        pieces = []
        for tx_id, features in genes[gene_id].items():
            by_locus = defaultdict(list)
            for f in features:
                by_locus[(f[0], f[6])].append(f)
            for (chrom, strand), feats in by_locus.items():
                start = min(int(f[3]) for f in feats)
                end = max(int(f[4]) for f in feats)
                pieces.append((chrom, strand, start, end, tx_id, feats))

        # Cluster pieces into loci: same sequence and strand, overlapping spans
        clusters = []
        for piece in sorted(pieces, key=lambda p: (p[0], p[1], p[2])):
            last = clusters[-1] if clusters else None
            if last and last['chrom'] == piece[0] and last['strand'] == piece[1] \
                    and piece[2] <= last['end']:
                last['pieces'].append(piece)
                last['end'] = max(last['end'], piece[3])
            else:
                clusters.append({'chrom': piece[0], 'strand': piece[1],
                                 'start': piece[2], 'end': piece[3],
                                 'pieces': [piece]})
        if len(clusters) < 2:
            continue

        clusters.sort(key=lambda c: (c['chrom'], c['start'], c['strand']))
        gene_template = gene_lines.get(gene_id)
        del genes[gene_id]
        new_gene_ids = []
        for k, cluster in enumerate(clusters):
            new_gene = gene_id if k == 0 else free_id(f"{gene_id}_{k + 1}", used_genes)
            new_gene_ids.append(new_gene)
            genes[new_gene] = {}
            for chrom, strand, _, _, tx_id, feats in cluster['pieces']:
                if k == 0 and tx_id not in genes[new_gene]:
                    new_tx = tx_id
                elif tx_id.startswith(gene_id + '.'):
                    new_tx = free_id(new_gene + tx_id[len(gene_id):], used_tx)
                else:
                    new_tx = free_id(f"{tx_id}_{k + 1}", used_tx)
                genes[new_gene][new_tx] = [
                    rename_ids(f, gene_id, new_gene, tx_id, new_tx) for f in feats]
                tx_to_gene[new_tx] = new_gene
                tx_template = tx_lines.get(tx_id)
                if tx_template:
                    tx_line = rename_ids(tx_template, gene_id, new_gene, tx_id, new_tx)
                    tx_line[0], tx_line[6] = chrom, strand
                    tx_lines[new_tx] = tx_line
            if gene_template:
                gene_line = rename_ids(gene_template, gene_id, new_gene)
                gene_line[0], gene_line[6] = cluster['chrom'], cluster['strand']
                gene_line[3], gene_line[4] = str(cluster['start']), str(cluster['end'])
                gene_lines[new_gene] = gene_line
        messages.append(f"SPLIT {gene_id}: {len(clusters)} loci -> "
                        f"{', '.join(new_gene_ids)}")

    return messages


def remove_duplicate_transcripts(genes, gene_lines, tx_lines):
    """Drop transcripts whose CDS chain repeats an earlier one (issue #71).

    Two transcripts are duplicates when they share sequence, strand and
    every CDS (start, end, frame), whether they belong to the same gene or
    not. The first one in genomic order is kept; a gene left without
    transcripts is dropped.

    Returns log messages.
    """
    messages = []
    seen = {}
    changed = set()

    def gene_key(gene_id):
        feats = [f for tx in genes[gene_id].values() for f in tx]
        return (feats[0][0], min(int(f[3]) for f in feats), gene_id) if feats \
            else ('', 0, gene_id)

    for gene_id in sorted(genes, key=gene_key):
        for tx_id in sorted(genes[gene_id]):
            cds = sorted((int(f[3]), int(f[4]), f[7])
                         for f in genes[gene_id][tx_id] if f[2] == 'CDS')
            if not cds:
                continue
            first = next(f for f in genes[gene_id][tx_id] if f[2] == 'CDS')
            key = (first[0], first[6], tuple(cds))
            if key in seen:
                messages.append(f"DUPLICATE {tx_id}: CDS identical to {seen[key]}, removed")
                del genes[gene_id][tx_id]
                changed.add(gene_id)
            else:
                seen[key] = tx_id
        if not genes[gene_id]:
            messages.append(f"DUPLICATE gene {gene_id}: no transcripts left, removed")
            del genes[gene_id]
            changed.discard(gene_id)

    for gene_id in changed:
        update_gene_boundaries(gene_id, genes, gene_lines, tx_lines)

    return messages


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

    # Give every locus its own gene ID before anything else
    all_messages = split_multilocus_ids(genes, gene_lines, tx_lines, tx_to_gene)
    split_count = len(all_messages)

    # Normalize each transcript
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

    # Drop repeated CDS chains (compared after stop codon trimming)
    dup_messages = remove_duplicate_transcripts(genes, gene_lines, tx_lines)
    all_messages.extend(dup_messages)
    dup_count = sum(1 for m in dup_messages if not m.startswith("DUPLICATE gene"))

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
        log.write(f"Gene IDs split (several loci): {split_count}\n")
        log.write(f"Duplicate transcripts removed (identical CDS): {dup_count}\n")
        log.write(f"Genes in output: {genes_after}\n")
        log.write(f"\nDetails:\n")
        for msg in all_messages:
            log.write(f"  {msg}\n")

    print(f"Done: {trimmed_count} trimmed, {len(discarded_genes)} discarded, "
          f"{non_atg_count} non-ATG starts, {split_count} gene IDs split, "
          f"{dup_count} duplicate transcripts removed, {genes_after} genes in output",
          file=sys.stderr)


if __name__ == "__main__":
    main()
