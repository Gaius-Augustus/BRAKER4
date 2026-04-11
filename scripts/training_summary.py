#!/usr/bin/env python3
"""
Generate AUGUSTUS training summary report and publication-quality plot.

Collects training gene counts at each pipeline stage and accuracy metrics
before and after AUGUSTUS parameter optimization. Produces:
1. A text summary file
2. A high-quality PDF/PNG bar+line plot

Usage:
    python3 training_summary.py -d output/sample/ -o output/sample/results/
"""

import argparse
import os
import re
import sys

def parse_accuracy_file(filepath):
    """Parse accuracy file, return dict with nu_sen, nu_sp, ex_sen, ex_sp, gene_sen, gene_sp, weighted."""
    acc = {}
    if not os.path.exists(filepath):
        return acc
    with open(filepath) as f:
        text = f.read()
    for label, key in [
        (r'Sensitivity:\s+([\d.]+)%', None),
        (r'Specificity:\s+([\d.]+)%', None),
    ]:
        pass  # use more specific patterns

    # Parse structured format
    section = None
    for line in text.splitlines():
        if 'Nucleotide level' in line:
            section = 'nu'
        elif 'Exon level' in line:
            section = 'ex'
        elif 'Gene level' in line:
            section = 'gene'
        elif 'Target accuracy' in line:
            m = re.search(r'([\d.]+)%', line)
            if m:
                acc['weighted'] = float(m.group(1))
        elif section:
            m = re.search(r'Sensitivity:\s+([\d.]+)%', line)
            if m:
                acc[f'{section}_sen'] = float(m.group(1))
            m = re.search(r'Specificity:\s+([\d.]+)%', line)
            if m:
                acc[f'{section}_sp'] = float(m.group(1))
    return acc


def parse_gene_count(filepath, pattern=None):
    """Extract a gene count from a file. Returns int or None."""
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        text = f.read().strip()
    if pattern:
        m = re.search(pattern, text)
        return int(m.group(1)) if m else None
    # Try "N -> M" format
    m = re.match(r'(\d+)\s*->\s*(\d+)', text)
    if m:
        return int(m.group(2))
    # Try just a number
    m = re.match(r'(\d+)', text)
    if m:
        return int(m.group(1))
    return None


def count_loci_in_gb(filepath):
    """Count LOCUS entries in GenBank file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        return sum(1 for line in f if line.startswith('LOCUS'))


def main():
    parser = argparse.ArgumentParser(description='Generate AUGUSTUS training summary')
    parser.add_argument('-d', '--workdir', required=True, help='Sample working directory (output/{sample}/)')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory for summary and plot')
    parser.add_argument('--no-plot', action='store_true', help='Skip plot generation')
    args = parser.parse_args()

    d = args.workdir
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    # Collect gene counts at each stage
    stages = []

    # Detect dual mode (two separate GeneMark-ETP runs that are merged)
    is_dual = (os.path.exists(os.path.join(d, 'GeneMark-ETP')) and
               os.path.exists(os.path.join(d, 'GeneMark-ETP-isoseq')))

    def count_unique_genes(gtf_path):
        """Count unique gene_id values in a GTF file.

        This is more robust than counting 'gene' feature lines because
        training.gtf files from GeneMark-ETP don't always include 'gene'
        features — they only have CDS/exon/intron lines with gene_id
        attributes.
        """
        if not os.path.exists(gtf_path):
            return 0
        gene_ids = set()
        with open(gtf_path) as f:
            for l in f:
                if l.startswith('#') or not l.strip():
                    continue
                m = re.search(r'gene_id "([^"]+)"', l)
                if m:
                    gene_ids.add(m.group(1))
        return len(gene_ids)

    if is_dual:
        # Dual mode: two parallel GeneMark-ETP runs, then merged
        n_sr_pred = count_unique_genes(os.path.join(d, 'GeneMark-ETP', 'genemark.gtf'))
        n_iso_pred = count_unique_genes(os.path.join(d, 'GeneMark-ETP-isoseq', 'genemark.gtf'))
        if n_sr_pred > 0:
            stages.append(('GeneMark-ETP (short-read)', n_sr_pred))
        if n_iso_pred > 0:
            stages.append(('GeneMark-ETP (IsoSeq)', n_iso_pred))

        # Training genes from each run
        n_sr = count_unique_genes(os.path.join(d, 'GeneMark-ETP', 'training.gtf'))
        n_iso = count_unique_genes(os.path.join(d, 'GeneMark-ETP-isoseq', 'training.gtf'))
        if n_sr > 0:
            stages.append(('Training genes (short-read)', n_sr))
        if n_iso > 0:
            stages.append(('Training genes (IsoSeq)', n_iso))

        # Merged + deduplicated training set
        merged_training = os.path.join(d, 'dual_etp_merged', 'training.gtf')
        if os.path.exists(merged_training):
            n_merged = count_unique_genes(merged_training)
            stages.append(('Training genes (merged)', n_merged))
    else:
        # Single-mode: one GeneMark run
        for gm_dir in ['genemark', 'GeneMark-ET', 'GeneMark-EP', 'GeneMark-ES', 'GeneMark-ETP']:
            gtf = os.path.join(d, gm_dir, 'genemark.gtf')
            if os.path.exists(gtf):
                n = count_unique_genes(gtf)
                stages.append(('GeneMark predictions', n))
                break

        # GeneMark filtered (good genes)
        for p in ['genemark/genemark.f.good.gtf', 'GeneMark-EP/genemark.f.good.gtf']:
            gtf = os.path.join(d, p)
            if os.path.exists(gtf):
                n = count_unique_genes(gtf)
                stages.append(('Filtered for training', n))
                break

        # ETP training genes (single-run)
        n = count_unique_genes(os.path.join(d, 'GeneMark-ETP', 'training.gtf'))
        if n > 0:
            stages.append(('ETP training genes', n))

    # After redundancy removal
    n = count_loci_in_gb(os.path.join(d, 'bonafide.f.gb'))
    if n is not None:
        stages.append(('After redundancy removal', n))

    # After etraining filtering
    n = count_loci_in_gb(os.path.join(d, 'bonafide.f.clean.gb'))
    if n is not None:
        stages.append(('After etraining filter', n))

    # After downsampling
    n = count_loci_in_gb(os.path.join(d, 'bonafide.f.clean.d.gb'))
    if n is not None:
        stages.append(('After downsampling', n))

    # Train/test split
    n_train = count_loci_in_gb(os.path.join(d, 'train.gb.train'))
    n_test = count_loci_in_gb(os.path.join(d, 'train.gb.test'))
    if n_train is not None:
        stages.append(('Training set', n_train))
    if n_test is not None:
        stages.append(('Test set (held out)', n_test))

    # Accuracy metrics
    acc_before = parse_accuracy_file(os.path.join(d, 'accuracy_after_training.txt'))
    acc_after = parse_accuracy_file(os.path.join(d, 'accuracy_after_optimize.txt'))

    # Final gene count
    braker_gtf = os.path.join(d, 'braker.gtf')
    if os.path.exists(braker_gtf):
        with open(braker_gtf) as f:
            n_final = sum(1 for line in f if '\tgene\t' in line)
    else:
        n_final = None

    # Write text summary
    summary_path = os.path.join(out, 'training_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("BRAKER AUGUSTUS Training Summary\n")
        f.write("=" * 50 + "\n\n")

        f.write("Training Gene Counts at Each Stage\n")
        f.write("-" * 50 + "\n")
        for label, count in stages:
            f.write(f"  {label:30s}  {count:>6d}\n")
        f.write("\n")

        if acc_before:
            f.write("Accuracy Before Optimization\n")
            f.write("-" * 50 + "\n")
            for level in ['nu', 'ex', 'gene']:
                name = {'nu': 'Nucleotide', 'ex': 'Exon', 'gene': 'Gene'}[level]
                sen = acc_before.get(f'{level}_sen', 0)
                sp = acc_before.get(f'{level}_sp', 0)
                f.write(f"  {name:12s}  Sn={sen:5.1f}%  Sp={sp:5.1f}%\n")
            f.write(f"  {'Weighted':12s}  {acc_before.get('weighted', 0):5.1f}%\n\n")

        if acc_after:
            f.write("Accuracy After Optimization\n")
            f.write("-" * 50 + "\n")
            for level in ['nu', 'ex', 'gene']:
                name = {'nu': 'Nucleotide', 'ex': 'Exon', 'gene': 'Gene'}[level]
                sen = acc_after.get(f'{level}_sen', 0)
                sp = acc_after.get(f'{level}_sp', 0)
                f.write(f"  {name:12s}  Sn={sen:5.1f}%  Sp={sp:5.1f}%\n")
            f.write(f"  {'Weighted':12s}  {acc_after.get('weighted', 0):5.1f}%\n")

            if acc_before and 'weighted' in acc_before and 'weighted' in acc_after:
                improvement = acc_after['weighted'] - acc_before['weighted']
                f.write(f"\n  Improvement: {improvement:+.1f}%\n")
            f.write("\n")

        if n_final is not None:
            f.write(f"Final prediction: {n_final} genes\n")

    print(f"Training summary written to {summary_path}")

    # Generate plot
    if args.no_plot:
        return

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
    except ImportError:
        print("WARNING: matplotlib not available, skipping plot", file=sys.stderr)
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [3, 2]})

    # Left panel: Training gene funnel
    if stages:
        labels = [s[0] for s in stages]
        counts = [s[1] for s in stages]
        colors = ['#4C72B0'] * len(stages)
        # Highlight by category
        for i, label in enumerate(labels):
            if 'Test' in label:
                colors[i] = '#DD8452'
            elif 'Training set' in label:
                colors[i] = '#55A868'
            elif 'short-read' in label:
                colors[i] = '#6A8CAE'  # lighter blue for short-read branch
            elif 'IsoSeq' in label:
                colors[i] = '#B08A7C'  # brown for IsoSeq branch
            elif 'merged' in label:
                colors[i] = '#8E7CC3'  # purple for merged result

        y_pos = range(len(labels))
        bars = ax1.barh(y_pos, counts, color=colors, edgecolor='white', linewidth=0.5)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(labels, fontsize=9)
        ax1.invert_yaxis()
        ax1.set_xlabel('Number of genes', fontsize=10)
        ax1.set_title('Training Gene Counts', fontsize=12, fontweight='bold')
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

        # Add count labels on bars
        for bar, count in zip(bars, counts):
            ax1.text(bar.get_width() + max(counts) * 0.02, bar.get_y() + bar.get_height() / 2,
                     str(count), va='center', fontsize=9)

    # Right panel: Accuracy comparison
    if acc_before or acc_after:
        metrics = ['Nuc Sn', 'Nuc Sp', 'Exon Sn', 'Exon Sp',
                   'Gene Sn', 'Gene Sp', 'Weighted']
        keys = ['nu_sen', 'nu_sp', 'ex_sen', 'ex_sp', 'gene_sen', 'gene_sp', 'weighted']

        before_vals = [acc_before.get(k, 0) for k in keys]
        after_vals = [acc_after.get(k, 0) for k in keys]

        x = range(len(metrics))
        width = 0.35

        if acc_before:
            ax2.bar([i - width / 2 for i in x], before_vals, width,
                    label='Before optimization', color='#C44E52', alpha=0.8)
        if acc_after:
            ax2.bar([i + width / 2 for i in x], after_vals, width,
                    label='After optimization', color='#55A868', alpha=0.8)

        ax2.set_ylabel('Accuracy (%)', fontsize=10)
        ax2.set_title('AUGUSTUS Training Accuracy', fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(metrics, fontsize=8, rotation=45, ha='right')
        ax2.set_ylim(0, 105)
        ax2.legend(fontsize=8, loc='lower right')
        ax2.axhline(y=50, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    plt.tight_layout()

    for fmt in ['pdf', 'png']:
        plot_path = os.path.join(out, f'training_summary.{fmt}')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Training plot saved to {os.path.join(out, 'training_summary.pdf')}")
    print(f"Training plot saved to {os.path.join(out, 'training_summary.png')}")


if __name__ == '__main__':
    main()
