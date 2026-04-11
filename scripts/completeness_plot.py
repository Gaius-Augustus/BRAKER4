#!/usr/bin/env python3
"""
Generate a combined BUSCO + compleasm completeness bar chart.

Shows 4 horizontal bars:
  - BUSCO genome
  - BUSCO proteome
  - compleasm genome
  - compleasm proteome

Each bar shows: Complete Single (S), Complete Duplicated (D), Fragmented (F), Missing (M)

Usage:
    python3 completeness_plot.py -d <results_qc_dir> -o <output.png>
"""

import argparse
import os
import re
import sys


def parse_busco_summary(path):
    """Parse BUSCO summary for genome and proteome scores."""
    results = {}
    if not os.path.exists(path):
        return results

    current = None
    with open(path) as f:
        for line in f:
            if "Genome" in line:
                current = "BUSCO genome"
            elif "Proteome" in line:
                current = "BUSCO proteome"
            elif current and "C:" in line:
                m = re.search(
                    r'C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%',
                    line
                )
                if m:
                    results[current] = {
                        'S': float(m.group(2)),
                        'D': float(m.group(3)),
                        'F': float(m.group(4)),
                        'M': float(m.group(5)),
                    }
                current = None
    return results


def parse_compleasm_summary(path, label):
    """Parse compleasm summary.txt."""
    if not os.path.exists(path):
        return None

    scores = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("S:"):
                scores['S'] = float(re.search(r'S:([0-9.]+)%', line).group(1))
            elif line.startswith("D:"):
                scores['D'] = float(re.search(r'D:([0-9.]+)%', line).group(1))
            elif line.startswith("F:"):
                scores['F'] = float(re.search(r'F:([0-9.]+)%', line).group(1))
            elif line.startswith("M:"):
                scores['M'] = float(re.search(r'M:([0-9.]+)%', line).group(1))

    if all(k in scores for k in ('S', 'D', 'F', 'M')):
        return {label: scores}
    return None


def generate_plot(data, output_path):
    """Generate horizontal stacked bar chart."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available, skipping plot", file=sys.stderr)
        return False

    labels = list(data.keys())
    S = [data[l]['S'] for l in labels]
    D = [data[l]['D'] for l in labels]
    F = [data[l]['F'] for l in labels]
    M = [data[l]['M'] for l in labels]

    fig, ax = plt.subplots(figsize=(10, max(2.5, len(labels) * 0.7 + 1)))

    y = np.arange(len(labels))
    bar_height = 0.5

    colors = {
        'S': '#2196F3',  # blue - single
        'D': '#64B5F6',  # light blue - duplicated
        'F': '#FFC107',  # amber - fragmented
        'M': '#F44336',  # red - missing
    }

    bars_s = ax.barh(y, S, bar_height, label='Complete (S)', color=colors['S'])
    bars_d = ax.barh(y, D, bar_height, left=S, label='Duplicated (D)', color=colors['D'])
    left_f = [s + d for s, d in zip(S, D)]
    bars_f = ax.barh(y, F, bar_height, left=left_f, label='Fragmented (F)', color=colors['F'])
    left_m = [l + f for l, f in zip(left_f, F)]
    bars_m = ax.barh(y, M, bar_height, left=left_m, label='Missing (M)', color=colors['M'])

    # Add percentage labels on bars
    for i in range(len(labels)):
        # Complete (S+D) label
        complete = S[i] + D[i]
        if complete > 8:
            ax.text(complete / 2, i, f"C:{complete:.1f}%",
                    ha='center', va='center', fontsize=9, fontweight='bold', color='white')

        # Missing label
        if M[i] > 8:
            ax.text(left_m[i] + M[i] / 2, i, f"M:{M[i]:.1f}%",
                    ha='center', va='center', fontsize=9, fontweight='bold', color='white')

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percentage of BUSCOs", fontsize=11)
    ax.legend(loc='lower right', fontsize=9, ncol=4)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate BUSCO/compleasm completeness plot")
    parser.add_argument("-d", "--qcdir", required=True, help="Quality control directory")
    parser.add_argument("-o", "--output", required=True, help="Output PNG path")
    args = parser.parse_args()

    data = {}

    # BUSCO
    busco_path = os.path.join(args.qcdir, "busco_summary.txt")
    busco = parse_busco_summary(busco_path)
    data.update(busco)

    # compleasm genome
    compleasm_genome = os.path.join(args.qcdir, "compleasm_genome_out", "summary.txt")
    result = parse_compleasm_summary(compleasm_genome, "compleasm genome")
    if result:
        data.update(result)

    # compleasm proteome
    compleasm_prot = os.path.join(args.qcdir, "compleasm_summary.txt")
    result = parse_compleasm_summary(compleasm_prot, "compleasm proteome")
    if result:
        data.update(result)

    if not data:
        print("No BUSCO/compleasm data found, skipping plot", file=sys.stderr)
        sys.exit(0)

    if generate_plot(data, args.output):
        print(f"Plot: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
