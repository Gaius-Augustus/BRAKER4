#!/usr/bin/env python3
"""Rescue multi-exon genes where AUGUSTUS and GeneMark agree but TSEBRA dropped them.

When both gene finders independently predict the same CDS intron chain at
a locus and TSEBRA still filtered the gene out, the agreement between two
independent methods is a strong confidence signal. This script adds those
genes back, provided every intron is supported by at least one hint.

Only multi-exon genes are rescued (single-exon agreement is vacuous).
Only genes NOT already in the BRAKER set are added (no modifications).

Benchmarked on A. thaliana ETP: +29 TP, +4 FP.
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_intron_hints(hints_path):
    """Parse hintsfile.gff and return a set of (chrom, start, end, strand) for intron hints.

    Gracefully handles files that contain no intron features (e.g. ES mode
    where the hints file is a dummy dependency).
    """
    intron_hints = set()
    try:
        with open(hints_path) as fh:
            for line in fh:
                if line.startswith('#') or line.strip() == '':
                    continue
                cols = line.split('\t')
                if len(cols) < 9:
                    continue
                if cols[2] != 'intron':
                    continue
                try:
                    chrom = cols[0]
                    start = int(cols[3])
                    end = int(cols[4])
                    strand = cols[6]
                    intron_hints.add((chrom, start, end, strand))
                except (ValueError, IndexError):
                    continue
    except OSError as e:
        print(f"[WARN] Could not read hints file {hints_path}: {e}", file=sys.stderr)
    return intron_hints


def parse_cds_intron_chains(gtf_path):
    """Parse a GTF and return per-transcript CDS intron chains.

    Returns:
        chains: dict mapping (chrom, strand, intron_chain_tuple) -> list of transcript_ids
        tx_cds: dict mapping transcript_id -> sorted list of (start, end) CDS tuples
        tx_info: dict mapping transcript_id -> {'chrom': str, 'strand': str}
    """
    tx_cds = defaultdict(list)
    tx_info = {}

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            if cols[2] != 'CDS':
                continue
            m = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
            if not m:
                continue
            tid = m.group(1)
            tx_cds[tid].append((int(cols[3]), int(cols[4])))
            tx_info[tid] = {'chrom': cols[0], 'strand': cols[6]}

    # Sort CDS and compute intron chains
    chains = defaultdict(list)
    for tid, cds_list in tx_cds.items():
        cds_list.sort()
        if len(cds_list) < 2:
            continue  # skip single-exon
        introns = []
        for i in range(len(cds_list) - 1):
            introns.append((cds_list[i][1] + 1, cds_list[i + 1][0] - 1))
        key = (tx_info[tid]['chrom'], tx_info[tid]['strand'], tuple(introns))
        chains[key].append(tid)

    return chains, tx_cds, tx_info


def extract_gene_blocks(gtf_path, target_tids):
    """Extract full gene blocks (all GTF lines) for given transcript IDs.

    Returns list of lines (strings) for each gene that contains at least
    one target transcript. Includes gene, transcript, CDS, exon, intron,
    start_codon, stop_codon, and all other feature lines.
    """
    target_tids = set(target_tids)

    # First pass: find gene_ids that contain target transcripts
    target_gids = set()
    for tid in target_tids:
        # gene_id is typically the part before the last dot
        # but let's parse it properly
        pass

    # Parse gene_id for each target transcript
    tid_to_gid = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            m_tid = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
            m_gid = re.search(r'gene_id\s+"([^"]+)"', cols[8])
            if m_tid and m_tid.group(1) in target_tids and m_gid:
                tid_to_gid[m_tid.group(1)] = m_gid.group(1)
                target_gids.add(m_gid.group(1))

    # For gene/transcript lines without standard attributes, check col[8] directly
    # AUGUSTUS format: gene line has just "g1" in col[8], transcript line has "g1.t1"

    # Second pass: collect all lines for target genes
    gene_lines = []
    in_target_gene = False
    current_gid = None

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue

            feature = cols[2]

            if feature == 'gene':
                # Gene line: col[8] is just the gene_id (AUGUSTUS format)
                gid = cols[8].strip().rstrip(';')
                # Also check for gene_id attribute format
                m_gid = re.search(r'gene_id\s+"([^"]+)"', cols[8])
                if m_gid:
                    gid = m_gid.group(1)
                if gid in target_gids:
                    in_target_gene = True
                    current_gid = gid
                else:
                    in_target_gene = False
                    current_gid = None

            if in_target_gene:
                gene_lines.append(line.rstrip('\n'))

    return gene_lines


def rename_gene_ids(lines, prefix="both_agree_"):
    """Rename gene and transcript IDs in GTF lines to avoid collisions.

    Prefixes all gene_id and transcript_id values. Also handles the bare
    AUGUSTUS-style gene/transcript lines where col[8] is just the ID.
    """
    renamed = []
    for line in lines:
        cols = line.split('\t')
        if len(cols) < 9:
            renamed.append(line)
            continue

        feature = cols[2]
        attr = cols[8]

        if feature == 'gene':
            # AUGUSTUS bare format: col[8] is just "g123"
            m = re.search(r'gene_id\s+"([^"]+)"', attr)
            if m:
                attr = attr.replace(f'"{m.group(1)}"', f'"{prefix}{m.group(1)}"')
            else:
                attr = prefix + attr
        elif feature == 'transcript':
            m_gid = re.search(r'gene_id\s+"([^"]+)"', attr)
            m_tid = re.search(r'transcript_id\s+"([^"]+)"', attr)
            if m_tid:
                attr = attr.replace(f'"{m_tid.group(1)}"', f'"{prefix}{m_tid.group(1)}"')
                if m_gid:
                    attr = attr.replace(f'"{m_gid.group(1)}"', f'"{prefix}{m_gid.group(1)}"')
            else:
                # Bare format: col[8] is "g123.t1"
                attr = prefix + attr
        else:
            # Standard feature lines with transcript_id and gene_id
            m_tid = re.search(r'transcript_id\s+"([^"]+)"', attr)
            m_gid = re.search(r'gene_id\s+"([^"]+)"', attr)
            if m_tid:
                attr = attr.replace(
                    f'transcript_id "{m_tid.group(1)}"',
                    f'transcript_id "{prefix}{m_tid.group(1)}"')
            if m_gid:
                attr = attr.replace(
                    f'gene_id "{m_gid.group(1)}"',
                    f'gene_id "{prefix}{m_gid.group(1)}"')

        cols[8] = attr
        renamed.append('\t'.join(cols))

    return renamed


def main():
    parser = argparse.ArgumentParser(
        description='Rescue multi-exon genes where AUGUSTUS and GeneMark agree but TSEBRA dropped.')
    parser.add_argument('--augustus', required=True,
                        help='AUGUSTUS hints GTF (augustus.hints.fixed.gtf)')
    parser.add_argument('--genemark', required=True,
                        help='GeneMark GTF')
    parser.add_argument('--braker', required=True,
                        help='BRAKER TSEBRA GTF (braker.tsebra.gtf)')
    parser.add_argument('--hints', required=True,
                        help='Merged hints file (hintsfile.gff)')
    parser.add_argument('--output', required=True,
                        help='Output GTF (braker + rescued genes)')
    args = parser.parse_args()

    # Step 1: Parse intron hints
    print("[INFO] Parsing intron hints...", file=sys.stderr)
    intron_hints = parse_intron_hints(args.hints)
    print(f"[INFO] Found {len(intron_hints)} unique intron hint coordinates", file=sys.stderr)

    # Step 2: Parse CDS intron chains from all three GTFs
    print("[INFO] Parsing AUGUSTUS GTF...", file=sys.stderr)
    aug_chains, aug_cds, aug_info = parse_cds_intron_chains(args.augustus)
    print(f"[INFO] AUGUSTUS: {len(aug_chains)} unique multi-exon intron chains", file=sys.stderr)

    print("[INFO] Parsing GeneMark GTF...", file=sys.stderr)
    gm_chains, _, _ = parse_cds_intron_chains(args.genemark)
    print(f"[INFO] GeneMark: {len(gm_chains)} unique multi-exon intron chains", file=sys.stderr)

    print("[INFO] Parsing BRAKER GTF...", file=sys.stderr)
    br_chains, _, _ = parse_cds_intron_chains(args.braker)
    print(f"[INFO] BRAKER: {len(br_chains)} unique multi-exon intron chains", file=sys.stderr)

    # Step 3: Find chains where AUGUSTUS and GeneMark agree but BRAKER lacks
    aug_keys = set(aug_chains.keys())
    gm_keys = set(gm_chains.keys())
    br_keys = set(br_chains.keys())

    both_agree = aug_keys & gm_keys        # both tools predicted this chain
    tsebra_dropped = both_agree - br_keys   # TSEBRA didn't keep it

    print(f"[INFO] Both-agree intron chains: {len(both_agree)}", file=sys.stderr)
    print(f"[INFO] Of those, TSEBRA dropped: {len(tsebra_dropped)}", file=sys.stderr)

    # Step 4: Filter by hint support -- every intron must have a hint
    rescue_tids = []
    n_no_hint = 0

    for key in tsebra_dropped:
        chrom, strand, intron_chain = key
        all_supported = True
        for intron_start, intron_end in intron_chain:
            if (chrom, intron_start, intron_end, strand) not in intron_hints:
                all_supported = False
                break
        if all_supported:
            # Pick the first AUGUSTUS transcript with this chain
            rescue_tids.append(aug_chains[key][0])
        else:
            n_no_hint += 1

    print(f"[INFO] Candidates passing hint filter: {len(rescue_tids)}", file=sys.stderr)
    print(f"[INFO] Candidates rejected (missing hint support): {n_no_hint}", file=sys.stderr)

    if not rescue_tids:
        print("[INFO] Nothing to rescue. Output will be identical to input BRAKER GTF.",
              file=sys.stderr)

    # Step 5: Extract full gene blocks from AUGUSTUS GTF and rename IDs
    if rescue_tids:
        print(f"[INFO] Extracting {len(rescue_tids)} gene blocks from AUGUSTUS GTF...",
              file=sys.stderr)
        gene_lines = extract_gene_blocks(args.augustus, rescue_tids)
        gene_lines = rename_gene_ids(gene_lines)
    else:
        gene_lines = []

    # Step 6: Write output = original BRAKER + rescued genes
    n_braker_lines = 0
    with open(args.braker) as fin, open(args.output, 'w') as fout:
        for line in fin:
            fout.write(line)
            n_braker_lines += 1
        if gene_lines:
            fout.write(f'# both-agree rescue: {len(rescue_tids)} genes added\n')
            for gline in gene_lines:
                fout.write(gline + '\n')

    print(f"[DONE] Wrote {args.output}: {n_braker_lines} original lines + "
          f"{len(gene_lines)} rescued lines ({len(rescue_tids)} genes)",
          file=sys.stderr)


if __name__ == '__main__':
    main()
