#!/usr/bin/env python3
"""
Generate comprehensive gene set statistics and publication-quality plots.

Analyzes braker.gtf to produce:
1. Summary statistics (genes, transcripts, exons, introns)
2. Transcripts-per-gene histogram
3. Single-exon vs multi-exon pie chart
4. Transcript length distributions (genomic span and CDS length)
5. Introns-per-gene histogram
6. Evidence support visualization (from gene_support.tsv)

All plots are saved as both PDF and PNG (300 DPI).

Usage:
    python3 gene_set_statistics.py -g braker.gtf -o results/quality_control/ \
        [-s gene_support.tsv] [-a braker.aa]
"""

import argparse
import os
import re
import sys
from collections import defaultdict, Counter


def parse_gtf(gtf_path):
    """Parse GTF file into gene/transcript/exon structure.

    Returns:
        genes: dict of gene_id -> {transcripts: {tx_id -> {exons: [(start,end,strand)], cds: [(s,e)]}}}
    """
    genes = defaultdict(lambda: {"transcripts": defaultdict(lambda: {"exons": [], "cds": []})})

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            feature = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]

            # Extract gene_id and transcript_id
            gm = re.search(r'gene_id "([^"]+)"', attrs)
            tm = re.search(r'transcript_id "([^"]+)"', attrs)

            if not gm:
                # Bare format: gene line has just gene_id in col9
                if feature == "gene":
                    gene_id = attrs.strip().rstrip(";")
                    genes[gene_id]  # ensure exists
                continue

            gene_id = gm.group(1)

            if feature == "gene":
                genes[gene_id]
            elif feature in ("transcript", "mRNA"):
                tx_id = tm.group(1) if tm else attrs.strip().rstrip(";")
                genes[gene_id]["transcripts"][tx_id]
            elif feature == "exon" and tm:
                tx_id = tm.group(1)
                genes[gene_id]["transcripts"][tx_id]["exons"].append((start, end, strand))
            elif feature == "CDS" and tm:
                tx_id = tm.group(1)
                genes[gene_id]["transcripts"][tx_id]["cds"].append((start, end))

    return dict(genes)


def compute_statistics(genes):
    """Compute summary statistics from parsed gene structure."""
    stats = {}

    n_genes = len(genes)
    n_transcripts = sum(len(g["transcripts"]) for g in genes.values())

    # Transcripts per gene
    tx_per_gene = [len(g["transcripts"]) for g in genes.values()]
    tx_per_gene_counts = Counter(tx_per_gene)

    # Exon counts per transcript
    exon_counts = []
    single_exon_tx = 0
    multi_exon_tx = 0
    cds_lengths = []  # sum of CDS exon lengths (excluding introns)
    genomic_spans = []  # from first exon start to last exon end (including introns)
    introns_per_gene = []

    for gene_id, gene in genes.items():
        gene_introns = 0
        for tx_id, tx in gene["transcripts"].items():
            exons = sorted(tx["exons"], key=lambda x: x[0])
            n_exons = len(exons) if exons else len(tx["cds"])

            if n_exons <= 1:
                single_exon_tx += 1
            else:
                multi_exon_tx += 1

            exon_counts.append(n_exons)

            # CDS length (sum of CDS exon lengths)
            cds = tx["cds"] if tx["cds"] else exons
            cds_len = sum(e - s + 1 for s, e in cds)
            if cds_len > 0:
                cds_lengths.append(cds_len)

            # Genomic span
            all_coords = exons if exons else [(s, e, ".") for s, e in tx["cds"]]
            if all_coords:
                span_start = min(c[0] for c in all_coords)
                span_end = max(c[1] for c in all_coords)
                genomic_spans.append(span_end - span_start + 1)

            # Introns for this transcript
            if len(exons) > 1:
                gene_introns += len(exons) - 1

        if gene["transcripts"]:
            # Average introns across transcripts for this gene
            n_tx = len(gene["transcripts"])
            introns_per_gene.append(gene_introns // n_tx if n_tx > 0 else 0)

    stats["n_genes"] = n_genes
    stats["n_transcripts"] = n_transcripts
    stats["tx_per_gene"] = tx_per_gene
    stats["tx_per_gene_counts"] = tx_per_gene_counts
    stats["single_exon_tx"] = single_exon_tx
    stats["multi_exon_tx"] = multi_exon_tx
    stats["exon_counts"] = exon_counts
    stats["cds_lengths"] = cds_lengths
    stats["genomic_spans"] = genomic_spans
    stats["introns_per_gene"] = introns_per_gene
    stats["mean_tx_per_gene"] = sum(tx_per_gene) / len(tx_per_gene) if tx_per_gene else 0
    stats["mean_exons"] = sum(exon_counts) / len(exon_counts) if exon_counts else 0
    stats["mean_cds_length"] = sum(cds_lengths) / len(cds_lengths) if cds_lengths else 0
    stats["mean_genomic_span"] = sum(genomic_spans) / len(genomic_spans) if genomic_spans else 0

    return stats


def parse_support_tsv(tsv_path):
    """Parse gene_support.tsv for evidence support visualization."""
    if not os.path.exists(tsv_path):
        return None

    data = []
    header = None
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            cols = line.split("\t")
            if len(cols) >= len(header):
                row = dict(zip(header, cols))
                data.append(row)
    return data if data else None


def generate_plots(stats, outdir, support_data=None):
    """Generate all publication-quality plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib/numpy not available, skipping plots", file=sys.stderr)
        return {}

    plots = {}
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    })

    # --- 1+2. Isoform distribution + exon structure (side by side) ---
    has_tx = bool(stats["tx_per_gene"])
    has_exon = stats["single_exon_tx"] + stats["multi_exon_tx"] > 0
    if has_tx or has_exon:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))

        if has_tx:
            max_tx = min(max(stats["tx_per_gene"]), 15)
            bins = range(1, max_tx + 2)
            counts, edges, patches = ax1.hist(stats["tx_per_gene"], bins=bins,
                                              color="#4C72B0", edgecolor="white",
                                              align="left", rwidth=0.8)
            total = len(stats["tx_per_gene"])
            for count, patch in zip(counts, patches):
                if count > 0:
                    pct = count / total * 100
                    ax1.text(patch.get_x() + patch.get_width() / 2, count + total * 0.01,
                             f"{int(count)}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=6.5)
            ax1.set_xlabel("Transcripts per gene")
            ax1.set_ylabel("Number of genes")
            ax1.set_title("Isoform Distribution")
            ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

        if has_exon:
            sizes = [stats["single_exon_tx"], stats["multi_exon_tx"]]
            labels = [f"Single-exon\n({stats['single_exon_tx']})",
                      f"Multi-exon\n({stats['multi_exon_tx']})"]
            colors = ["#DD8452", "#4C72B0"]
            ax2.pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%",
                    startangle=90, textprops={"fontsize": 9})
            ax2.set_title(f"Exon Structure ({sum(sizes)} transcripts)")

        plt.tight_layout()
        path = os.path.join(outdir, "isoform_and_exon_structure.png")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close()
        plots["isoform_and_exon_structure"] = path

    # --- 3. Transcript length distributions (dual histogram) ---
    if stats["cds_lengths"] and stats["genomic_spans"]:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

        # CDS length (excluding introns)
        ax1.hist(stats["cds_lengths"], bins=50, color="#55A868", edgecolor="white", alpha=0.8)
        ax1.set_xlabel("CDS length (bp)")
        ax1.set_ylabel("Number of transcripts")
        ax1.set_title("CDS Length Distribution")
        median_cds = sorted(stats["cds_lengths"])[len(stats["cds_lengths"]) // 2]
        ax1.axvline(median_cds, color="red", linestyle="--", linewidth=1,
                    label=f"Median: {median_cds:,} bp")
        ax1.legend()

        # Genomic span (including introns)
        ax2.hist(stats["genomic_spans"], bins=50, color="#4C72B0", edgecolor="white", alpha=0.8)
        ax2.set_xlabel("Genomic span (bp)")
        ax2.set_ylabel("Number of transcripts")
        ax2.set_title("Genomic Span Distribution (incl. introns)")
        median_span = sorted(stats["genomic_spans"])[len(stats["genomic_spans"]) // 2]
        ax2.axvline(median_span, color="red", linestyle="--", linewidth=1,
                    label=f"Median: {median_span:,} bp")
        ax2.legend()

        plt.tight_layout()
        path = os.path.join(outdir, "transcript_lengths.png")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close()
        plots["transcript_lengths"] = path

    # --- 4. Introns per gene histogram (compact) ---
    if stats["introns_per_gene"]:
        fig, ax = plt.subplots(figsize=(5, 3))
        max_introns = min(max(stats["introns_per_gene"]) + 1, 30)
        ax.hist(stats["introns_per_gene"], bins=range(0, max_introns + 1),
                color="#C44E52", edgecolor="white", align="left", rwidth=0.8)
        ax.set_xlabel("Introns per gene")
        ax.set_ylabel("Number of genes")
        ax.set_title("Intron Count Distribution")
        median_introns = sorted(stats["introns_per_gene"])[len(stats["introns_per_gene"]) // 2]
        ax.axvline(median_introns, color="black", linestyle="--", linewidth=1,
                   label=f"Median: {median_introns}")
        ax.legend()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tight_layout()
        path = os.path.join(outdir, "introns_per_gene.png")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close()
        plots["introns_per_gene"] = path

    # --- 5. Evidence support visualization ---
    if support_data:
        # Parse support columns — look for intron/exon support counts
        try:
            # Typical columns: transcript_id, n_introns, n_supported_introns, ...
            col_names = list(support_data[0].keys())

            # Find support-related columns (handles both "support" and "sup" naming)
            intron_cols = [c for c in col_names if "intron" in c.lower() and ("support" in c.lower() or "sup" in c.lower())]
            exon_cols = [c for c in col_names if "exon" in c.lower() and ("support" in c.lower() or "sup" in c.lower())]
            # Find the total introns column (num_introns or n_introns)
            n_introns_col = next((c for c in col_names if c in ("n_introns", "num_introns")), None)

            if intron_cols or exon_cols:
                fig, axes = plt.subplots(1, 2, figsize=(12, 4))

                # Intron support fraction per transcript
                # Prefer "any" support column (combined RNA-Seq + protein)
                intron_sup_col = next((c for c in intron_cols if "any" in c.lower()), intron_cols[0] if intron_cols else None)
                if intron_sup_col and n_introns_col:
                    n_introns = [int(r.get(n_introns_col, 0)) for r in support_data]
                    n_supported = [int(r.get(intron_sup_col, 0)) for r in support_data]
                    fractions = [s / n if n > 0 else 0 for s, n in zip(n_supported, n_introns)]
                    axes[0].hist(fractions, bins=20, color="#4C72B0", edgecolor="white")
                    axes[0].set_xlabel("Fraction of introns supported")
                    axes[0].set_ylabel("Number of transcripts")
                    axes[0].set_title("Intron Evidence Support")
                    axes[0].set_xlim(-0.05, 1.05)

                    # Categorize by support level
                    full_support = sum(1 for f in fractions if f >= 1.0)
                    partial = sum(1 for f in fractions if 0 < f < 1.0)
                    none = sum(1 for f in fractions if f == 0)
                    total = len(fractions)

                    sizes = [full_support, partial, none]
                    labels = [f"Full support\n({full_support:,})",
                              f"Partial\n({partial:,})",
                              f"No support\n({none:,})"]
                    colors = ["#55A868", "#DD8452", "#C44E52"]
                    axes[1].pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%",
                               startangle=90)
                    axes[1].set_title(f"Evidence Support ({total:,} transcripts)")

                    plt.tight_layout()
                    path = os.path.join(outdir, "evidence_support.png")
                    fig.savefig(path, dpi=300, bbox_inches="tight")
                    fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
                    plt.close()
                    plots["evidence_support"] = path
        except (KeyError, ValueError, IndexError):
            pass  # Skip if support data format doesn't match

    return plots


def write_summary(stats, outdir):
    """Write gene set statistics summary text file."""
    path = os.path.join(outdir, "gene_set_statistics.txt")
    with open(path, "w") as f:
        f.write("BRAKER4 Gene Set Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Genes:                  {stats['n_genes']:>8,}\n")
        f.write(f"Transcripts:            {stats['n_transcripts']:>8,}\n")
        f.write(f"Mean transcripts/gene:  {stats['mean_tx_per_gene']:>8.2f}\n")
        f.write(f"Mean exons/transcript:  {stats['mean_exons']:>8.2f}\n")
        f.write(f"\n")
        f.write(f"Single-exon transcripts:{stats['single_exon_tx']:>8,} ({stats['single_exon_tx']/(stats['single_exon_tx']+stats['multi_exon_tx'])*100:.1f}%)\n")
        f.write(f"Multi-exon transcripts: {stats['multi_exon_tx']:>8,} ({stats['multi_exon_tx']/(stats['single_exon_tx']+stats['multi_exon_tx'])*100:.1f}%)\n")
        f.write(f"\n")
        f.write(f"Mean CDS length:        {stats['mean_cds_length']:>8,.0f} bp\n")
        f.write(f"Mean genomic span:      {stats['mean_genomic_span']:>8,.0f} bp\n")
        if stats["cds_lengths"]:
            median_cds = sorted(stats["cds_lengths"])[len(stats["cds_lengths"]) // 2]
            f.write(f"Median CDS length:      {median_cds:>8,} bp\n")
        if stats["genomic_spans"]:
            median_span = sorted(stats["genomic_spans"])[len(stats["genomic_spans"]) // 2]
            f.write(f"Median genomic span:    {median_span:>8,} bp\n")

    return path


def main():
    parser = argparse.ArgumentParser(description="Generate gene set statistics and plots")
    parser.add_argument("-g", "--gtf", required=True, help="BRAKER GTF file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-s", "--support", default=None, help="gene_support.tsv file")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Parsing {args.gtf}...")
    genes = parse_gtf(args.gtf)
    stats = compute_statistics(genes)

    print(f"  {stats['n_genes']} genes, {stats['n_transcripts']} transcripts")

    summary_path = write_summary(stats, args.outdir)
    print(f"Summary: {summary_path}")

    support_data = parse_support_tsv(args.support) if args.support else None

    plots = generate_plots(stats, args.outdir, support_data)
    for name, path in plots.items():
        print(f"Plot: {path}")


if __name__ == "__main__":
    main()
