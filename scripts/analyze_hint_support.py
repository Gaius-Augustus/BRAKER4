#!/usr/bin/env python3
"""
Analyze hint support for transcripts in a GTF file.

This script analyzes how well transcripts are supported by hints (evidence) by examining:
- Multi-exon transcripts: full support (all introns covered), partial support (some introns covered), no support
- Single-exon transcripts: CDS/CDSpart hint coverage vs no coverage
- Intron hint usage: total hints, incorporated hints, unincorporated hints with high multiplicity (mult > 10)

High-multiplicity unincorporated hints indicate strong evidence that wasn't used in gene models,
which may suggest missing or incomplete gene predictions.

Author: Katharina J. Hoff
Date: 2025-10-12
"""

import argparse
import sys
from collections import defaultdict


def parse_gtf_line(line):
    """Parse a GTF line and return relevant information."""
    if line.startswith('#') or not line.strip():
        return None

    fields = line.strip().split('\t')
    if len(fields) < 9:
        return None

    seqname, source, feature, start, end, score, strand, frame, attributes = fields

    # Parse attributes
    attr_dict = {}
    for attr in attributes.split(';'):
        attr = attr.strip()
        if not attr:
            continue
        if '"' in attr:
            key, value = attr.split(' ', 1)
            value = value.strip('"')
            attr_dict[key] = value

    return {
        'seqname': seqname,
        'source': source,
        'feature': feature,
        'start': int(start),
        'end': int(end),
        'score': score,
        'strand': strand,
        'frame': frame,
        'attributes': attr_dict,
        'gene_id': attr_dict.get('gene_id'),
        'transcript_id': attr_dict.get('transcript_id')
    }


def parse_gff_line(line):
    """Parse a GFF/hints line and return relevant information."""
    if line.startswith('#') or not line.strip():
        return None

    fields = line.strip().split('\t')
    if len(fields) < 9:
        return None

    seqname, source, feature, start, end, score, strand, frame, attributes = fields

    # Parse attributes to extract mult value
    mult = 1
    for attr in attributes.split(';'):
        attr = attr.strip()
        if attr.startswith('mult='):
            try:
                mult = int(attr.split('=')[1])
            except (ValueError, IndexError):
                mult = 1

    return {
        'seqname': seqname,
        'source': source,
        'feature': feature,
        'start': int(start),
        'end': int(end),
        'score': score,
        'strand': strand,
        'frame': frame,
        'attributes': attributes,
        'mult': mult
    }


def load_transcripts_and_introns(gtf_file):
    """
    Load transcripts and their introns from GTF file.

    Returns:
        transcripts: dict of transcript_id -> {
            'exon_count': int,
            'introns': list of (seqname, start, end, strand),
            'seqname': str,
            'strand': str,
            'gene_id': str
        }
    """
    transcripts = defaultdict(lambda: {
        'exon_count': 0,
        'introns': [],
        'exons': [],
        'seqname': None,
        'strand': None,
        'gene_id': None
    })

    print(f"[INFO] Loading transcripts and introns from {gtf_file}...", file=sys.stderr)

    with open(gtf_file, 'r') as f:
        for line in f:
            data = parse_gtf_line(line)
            if not data or not data['transcript_id']:
                continue

            tx_id = data['transcript_id']
            gene_id = data['gene_id']

            # Store basic transcript info
            if transcripts[tx_id]['seqname'] is None:
                transcripts[tx_id]['seqname'] = data['seqname']
                transcripts[tx_id]['strand'] = data['strand']
                transcripts[tx_id]['gene_id'] = gene_id

            # Count exons and store exon positions
            if data['feature'] == 'exon':
                transcripts[tx_id]['exon_count'] += 1
                transcripts[tx_id]['exons'].append({
                    'start': data['start'],
                    'end': data['end']
                })

    # Calculate introns from exons for each transcript
    for tx_id, tx_data in transcripts.items():
        exons = sorted(tx_data['exons'], key=lambda x: x['start'])

        if len(exons) < 2:
            continue

        # Introns are between consecutive exons
        for i in range(len(exons) - 1):
            intron_start = exons[i]['end'] + 1
            intron_end = exons[i + 1]['start'] - 1

            intron = (tx_data['seqname'], intron_start, intron_end, tx_data['strand'])
            transcripts[tx_id]['introns'].append(intron)

    # Convert defaultdict to regular dict and remove temporary exon list
    transcripts = dict(transcripts)
    for tx_id in transcripts:
        del transcripts[tx_id]['exons']

    print(f"[INFO] Loaded {len(transcripts)} transcripts", file=sys.stderr)

    return transcripts


def load_hints(hints_file):
    """
    Load hints from GFF file.

    Returns:
        intron_hints: set of (seqname, start, end, strand) for intron hints
        intron_hints_with_mult: dict of (seqname, start, end, strand) -> mult value
        cds_hints: dict of seqname -> list of (start, end, strand) for CDS/CDSpart hints
    """
    intron_hints = set()
    intron_hints_with_mult = {}
    cds_hints = defaultdict(list)

    print(f"[INFO] Loading hints from {hints_file}...", file=sys.stderr)

    with open(hints_file, 'r') as f:
        for line in f:
            data = parse_gff_line(line)
            if not data:
                continue

            if data['feature'] == 'intron':
                hint = (data['seqname'], data['start'], data['end'], data['strand'])
                intron_hints.add(hint)
                intron_hints_with_mult[hint] = data['mult']
            elif data['feature'] in ('CDS', 'CDSpart'):
                cds_hints[data['seqname']].append({
                    'start': data['start'],
                    'end': data['end'],
                    'strand': data['strand']
                })

    print(f"[INFO] Loaded {len(intron_hints)} intron hints and {sum(len(v) for v in cds_hints.values())} CDS hints", file=sys.stderr)

    return intron_hints, intron_hints_with_mult, cds_hints


def check_intron_support(intron, intron_hints):
    """Check if an intron is supported by a hint."""
    return intron in intron_hints


def check_cds_overlap(tx_seqname, tx_strand, cds_hints):
    """Check if there's any CDS hint overlap for a transcript's location."""
    if tx_seqname not in cds_hints:
        return False

    # For simplicity, just check if there's any CDS hint on the same sequence and strand
    for hint in cds_hints[tx_seqname]:
        if hint['strand'] == tx_strand:
            return True

    return False


def analyze_hint_support(transcripts, intron_hints, intron_hints_with_mult, cds_hints):
    """
    Analyze hint support for all transcripts.

    Returns:
        stats: dict with various statistics
    """
    stats = {
        'multi_exon_full_support': 0,
        'multi_exon_partial_support': 0,
        'multi_exon_no_support': 0,
        'single_exon_with_cds_hint': 0,
        'single_exon_no_cds_hint': 0,
        'total_multi_exon': 0,
        'total_single_exon': 0,
        'total_transcripts': len(transcripts),
        'total_genes': len(set(tx['gene_id'] for tx in transcripts.values() if tx['gene_id']))
    }

    # Track which intron hints were incorporated into transcripts
    incorporated_hints = set()

    print(f"[INFO] Analyzing hint support...", file=sys.stderr)

    for tx_id, tx_data in transcripts.items():
        exon_count = tx_data['exon_count']
        introns = tx_data['introns']

        if exon_count == 1:
            # Single-exon transcript
            stats['total_single_exon'] += 1

            # Check for CDS hint support
            if check_cds_overlap(tx_data['seqname'], tx_data['strand'], cds_hints):
                stats['single_exon_with_cds_hint'] += 1
            else:
                stats['single_exon_no_cds_hint'] += 1

        elif exon_count > 1:
            # Multi-exon transcript
            stats['total_multi_exon'] += 1

            if len(introns) == 0:
                # No introns identified (shouldn't happen if exon_count > 1, but handle it)
                stats['multi_exon_no_support'] += 1
                continue

            # Check how many introns have hint support and track incorporated hints
            supported_introns = 0
            for intron in introns:
                if check_intron_support(intron, intron_hints):
                    supported_introns += 1
                    incorporated_hints.add(intron)

            if supported_introns == len(introns):
                stats['multi_exon_full_support'] += 1
            elif supported_introns > 0:
                stats['multi_exon_partial_support'] += 1
            else:
                stats['multi_exon_no_support'] += 1

    # Count unincorporated intron hints with high multiplicity
    unincorporated_mult_gt_10 = 0
    unincorporated_mult_gt_100 = 0
    for hint, mult in intron_hints_with_mult.items():
        if hint not in incorporated_hints:
            if mult > 10:
                unincorporated_mult_gt_10 += 1
            if mult > 100:
                unincorporated_mult_gt_100 += 1

    stats['unincorporated_hints_mult_gt_10'] = unincorporated_mult_gt_10
    stats['unincorporated_hints_mult_gt_100'] = unincorporated_mult_gt_100
    stats['total_intron_hints'] = len(intron_hints)
    stats['incorporated_hints'] = len(incorporated_hints)
    stats['unincorporated_hints'] = len(intron_hints) - len(incorporated_hints)

    return stats


def print_stats(stats, output_file=None):
    """Print statistics to stdout or file."""
    output = []

    output.append("=== Hint Support Analysis ===\n")
    output.append(f"Total genes: {stats['total_genes']}")
    output.append(f"Total transcripts: {stats['total_transcripts']}")
    output.append("")

    output.append("Multi-exon transcripts:")
    output.append(f"  Total: {stats['total_multi_exon']}")
    output.append(f"  Full support (all introns): {stats['multi_exon_full_support']} ({stats['multi_exon_full_support'] / max(1, stats['total_multi_exon']) * 100:.1f}%)")
    output.append(f"  Partial support (some introns): {stats['multi_exon_partial_support']} ({stats['multi_exon_partial_support'] / max(1, stats['total_multi_exon']) * 100:.1f}%)")
    output.append(f"  No support: {stats['multi_exon_no_support']} ({stats['multi_exon_no_support'] / max(1, stats['total_multi_exon']) * 100:.1f}%)")

    output.append("")
    output.append("Single-exon transcripts:")
    output.append(f"  Total: {stats['total_single_exon']}")
    output.append(f"  With CDS/CDSpart hint: {stats['single_exon_with_cds_hint']} ({stats['single_exon_with_cds_hint'] / max(1, stats['total_single_exon']) * 100:.1f}%)")
    output.append(f"  No CDS hint: {stats['single_exon_no_cds_hint']} ({stats['single_exon_no_cds_hint'] / max(1, stats['total_single_exon']) * 100:.1f}%)")

    output.append("")
    output.append("Intron hint usage:")
    output.append(f"  Total intron hints: {stats['total_intron_hints']}")
    output.append(f"  Incorporated into transcripts: {stats['incorporated_hints']} ({stats['incorporated_hints'] / max(1, stats['total_intron_hints']) * 100:.1f}%)")
    output.append(f"  Unincorporated: {stats['unincorporated_hints']} ({stats['unincorporated_hints'] / max(1, stats['total_intron_hints']) * 100:.1f}%)")
    output.append(f"    - Unincorporated with mult > 10: {stats['unincorporated_hints_mult_gt_10']} ({stats['unincorporated_hints_mult_gt_10'] / max(1, stats['unincorporated_hints']) * 100:.1f}% of unincorporated)")
    output.append(f"    - Unincorporated with mult > 100: {stats['unincorporated_hints_mult_gt_100']} ({stats['unincorporated_hints_mult_gt_100'] / max(1, stats['unincorporated_hints']) * 100:.1f}% of unincorporated)")

    output_text = '\n'.join(output) + '\n'

    if output_file:
        with open(output_file, 'w') as f:
            f.write(output_text)
        print(f"[INFO] Statistics written to {output_file}", file=sys.stderr)
    else:
        print(output_text)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze hint support for transcripts in a GTF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --gtf braker.gtf --hints hintsfile.gff
  %(prog)s --gtf braker.gtf --hints hintsfile.gff --output hint_stats.txt

Note: This script analyzes hint support per transcript, not per gene.
A gene may have multiple transcripts (isoforms) with different hint support.
        """
    )

    parser.add_argument('--gtf', required=True,
                        help='Input GTF file with gene predictions')
    parser.add_argument('--hints', required=True,
                        help='Input GFF file with hints (evidence)')
    parser.add_argument('--output', '-o',
                        help='Output file for statistics (default: stdout)')

    args = parser.parse_args()

    # Load data
    transcripts = load_transcripts_and_introns(args.gtf)
    intron_hints, intron_hints_with_mult, cds_hints = load_hints(args.hints)

    # Analyze
    stats = analyze_hint_support(transcripts, intron_hints, intron_hints_with_mult, cds_hints)

    # Print results
    print_stats(stats, args.output)

    print("[INFO] Analysis complete!", file=sys.stderr)


if __name__ == '__main__':
    main()
