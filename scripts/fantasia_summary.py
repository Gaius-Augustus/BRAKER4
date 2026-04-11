#!/usr/bin/env python3
"""
Summarise a FANTASIA-Lite results.csv into a text report and a GO-namespace plot.

FANTASIA-Lite emits one row per (protein, GO term) pair with the columns:
    query_accession, go_id, reliability_index, distance, go_description, category

`category` is the GO namespace (biological_process, molecular_function,
cellular_component). `reliability_index` is FANTASIA-Lite's confidence score
(higher is better; the per-team default cutoff is 0.5).

Outputs:
    - <out_dir>/fantasia_summary.txt          (plain-text key/value summary)
    - <out_dir>/fantasia_go_categories.png    (pie chart of broad functional
                                              categories; each protein
                                              contributes once per matching
                                              category)
    - <out_dir>/fantasia_go_terms.tsv         (flat per-(transcript, GO term)
                                              table with human-readable GO
                                              names; intended as the primary
                                              user-facing functional table)
"""

import argparse
import csv
import os
import re as _re
import sys
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


NAMESPACE_LABELS = {
    "biological_process": "Biological process",
    "molecular_function": "Molecular function",
    "cellular_component": "Cellular component",
    # FANTASIA-Lite emits single-letter namespace codes too:
    "P": "Biological process",
    "F": "Molecular function",
    "C": "Cellular component",
}


# ----------------------------------------------------------------------------
# Functional categories
# ----------------------------------------------------------------------------
# Broad biologically-meaningful buckets adapted from the EukAssembly-Bin
# functional_profile.py categorisation, trimmed to themes that apply to a
# single eukaryotic genome (metagenomics-specific cycles like nitrogen and
# sulfur cycling are dropped, eukaryote-specific themes like cytoskeleton
# and cell cycle are added). Each category is a set of "seed" GO IDs plus
# a list of regex patterns matched against the GO term name (FANTASIA-Lite
# provides this inline as `go_description`).
#
# Each protein is counted ONCE per category it matches (multi-category
# membership allowed); the pie chart therefore shows fraction of category
# memberships, not fraction of proteins.
FUNCTIONAL_CATEGORIES = [
    ("Energy / photosynthesis / respiration", {
        "go_ids": {"GO:0006091", "GO:0015979", "GO:0022900", "GO:0006119",
                   "GO:0009060", "GO:0042773", "GO:0009521", "GO:0009522",
                   "GO:0009523", "GO:0030096", "GO:0006099", "GO:0006101",
                   "GO:0006102", "GO:0006103"},
        "keywords": [r"electron.?transport", r"photosynth", r"respirat",
                     r"cytochrom", r"\bNADH\b", r"ATP.?synth",
                     r"chlorophyll", r"photosystem", r"thylakoid",
                     r"tricarboxylic", r"\bTCA\b", r"citric.acid.cycle",
                     r"isocitrate", r"\bcitrate\b", r"succinate",
                     r"fumarate", r"\bmalate\b", r"oxoglutarate",
                     r"oxidative.phosphoryl"],
    }),
    ("Carbohydrate metabolism", {
        "go_ids": {"GO:0005975", "GO:0016798", "GO:0004553", "GO:0016799",
                   "GO:0000272", "GO:0030245", "GO:0052689", "GO:0016052",
                   "GO:0005976", "GO:0016757", "GO:0030246", "GO:0030247"},
        "keywords": [r"glycos(yl|id)", r"carbohydrate", r"polysaccharide",
                     r"cellulos", r"chitin", r"starch", r"sucros", r"mannos"],
    }),
    ("Lipid metabolism", {
        "go_ids": {"GO:0006629", "GO:0008610", "GO:0006631", "GO:0006633",
                   "GO:0006635", "GO:0030148", "GO:0006643", "GO:0006644"},
        "keywords": [r"\blipid", r"fatty.?acid", r"phospholip", r"sphingolip",
                     r"sterol", r"triglycer", r"cholester"],
    }),
    ("Amino acid metabolism", {
        "go_ids": {"GO:0006520", "GO:0008652", "GO:0009063", "GO:0006518"},
        "keywords": [r"amino.?acid", r"aminotransferase",
                     r"glutamate", r"arginine", r"lysine", r"methionine"],
    }),
    ("Nucleotide metabolism", {
        "go_ids": {"GO:0009117", "GO:0006139", "GO:0009165"},
        "keywords": [r"nucleotid", r"\bpurin", r"pyrimidin", r"ribonucleotid"],
    }),
    ("DNA replication / repair / chromatin", {
        "go_ids": {"GO:0006260", "GO:0006281", "GO:0006289", "GO:0006298",
                   "GO:0003677", "GO:0003887", "GO:0006334", "GO:0000785"},
        "keywords": [r"DNA.replication", r"DNA.polymerase", r"DNA.repair",
                     r"helicase", r"topoisomerase", r"recombinase",
                     r"mismatch.repair", r"chromatin", r"nucleosome",
                     r"histone"],
    }),
    ("Transcription / RNA processing", {
        "go_ids": {"GO:0006351", "GO:0006355", "GO:0003700", "GO:0006366",
                   "GO:0006397", "GO:0008380", "GO:0003723", "GO:0003729"},
        "keywords": [r"transcription", r"RNA.polymerase",
                     r"transcription.?factor", r"RNA.splicing", r"\bspliceos",
                     r"mRNA.processing", r"RNA.binding", r"mRNA.binding",
                     r"rRNA.processing", r"snoRNA"],
    }),
    ("Translation / ribosome", {
        "go_ids": {"GO:0006412", "GO:0003735", "GO:0006414", "GO:0006418",
                   "GO:0005840"},
        "keywords": [r"ribosomal", r"\btranslation", r"\btRNA\b",
                     r"aminoacyl.tRNA", r"elongation.factor"],
    }),
    ("Proteolysis / protein folding", {
        "go_ids": {"GO:0006508", "GO:0008233", "GO:0004252", "GO:0004175",
                   "GO:0006457", "GO:0000502"},
        "keywords": [r"\bprotease", r"peptidase", r"proteasome",
                     r"ubiquitin", r"chaperon", r"protein.folding"],
    }),
    ("Transport", {
        "go_ids": {"GO:0006810", "GO:0055085", "GO:0015031", "GO:0015103",
                   "GO:0005215", "GO:0008643", "GO:0042626"},
        "keywords": [r"\btransport", r"\btransporter", r"ABC.?transport",
                     r"\bMFS\b", r"porin", r"\bchannel\b", r"permease",
                     r"symporter", r"antiporter"],
    }),
    ("Signal transduction", {
        "go_ids": {"GO:0007165", "GO:0023052", "GO:0000160", "GO:0004672",
                   "GO:0035556"},
        "keywords": [r"\bkinase", r"signal.transduc", r"\bGTPase",
                     r"phosphoryl(at|as)", r"phosphatase"],
    }),
    ("Stress / redox", {
        "go_ids": {"GO:0006950", "GO:0006979", "GO:0009408", "GO:0006970",
                   "GO:0004601"},
        "keywords": [r"stress.response", r"oxidative.stress", r"heat.shock",
                     r"catalase", r"superoxide", r"peroxidase",
                     r"thioredoxin", r"glutathione"],
    }),
    ("Cytoskeleton / cell structure", {
        "go_ids": {"GO:0005856", "GO:0005874", "GO:0005884", "GO:0005819",
                   "GO:0007010", "GO:0030863"},
        "keywords": [r"cytoskelet", r"microtubule", r"\bactin\b", r"tubulin",
                     r"intermediate.filament", r"cortical.cytoskel"],
    }),
    ("Cilium / flagellum / motility", {
        "go_ids": {"GO:0005929", "GO:0036157", "GO:0031514", "GO:0097014",
                   "GO:0035082", "GO:0001539", "GO:0003341"},
        "keywords": [r"\bcilium\b", r"\bcilia\b", r"\bciliary",
                     r"\bflagell", r"axoneme", r"intraflagellar",
                     r"basal.body", r"motile.cilium", r"dynein"],
    }),
    ("Cell wall / membrane", {
        "go_ids": {"GO:0071555", "GO:0009273", "GO:0016020", "GO:0042546",
                   "GO:0005886", "GO:0005618"},
        "keywords": [r"cell.?wall", r"plasma.membrane", r"lipopolysaccharide",
                     r"murein", r"peptidoglycan"],
    }),
    ("Cell cycle / division", {
        "go_ids": {"GO:0007049", "GO:0051301", "GO:0000278", "GO:0007059",
                   "GO:0000910"},
        "keywords": [r"cell.cycle", r"\bmitosis\b", r"\bmeiosis\b",
                     r"cytokinesis", r"chromosome.segregation",
                     r"spindle"],
    }),
]
_OTHER_LABEL = "Other / unclassified"
_LOC_ONLY_LABEL = "Subcellular localization only"

# Subcellular-localization GO term names. Used as a fallback bucket for
# proteins whose only annotations are "where they live" rather than "what
# they do" -- this is informative ("we know it's nuclear, just not what it
# does there") and keeps the residual Other / unclassified slice small.
_LOCALIZATION_KW = _re.compile(
    r"\b(?:cytoplasm|cytosol|nucleus|nucle(?:olus|oplasm|ar.lumen|ar.envelope|ar.pore)"
    r"|mitochondri[oa]|chloroplast|plastid|endoplasmic.reticulum|golgi"
    r"|peroxisome|lysosome|vacuole|vesicle|endosome|glyoxysome|glycosome"
    r"|kinetoplast|apicoplast|extracellular.space|extracellular.region"
    r"|integral.component.of.membrane|membrane.raft|cell.cortex)\b",
    _re.IGNORECASE,
)

# Distinct, mostly-colorblind-friendly palette for the pie wedges. Order is
# stable so re-runs produce identical figures across samples.
_CATEGORY_COLORS = [
    "#3a7d44", "#e07b00", "#2e7da8", "#c46a1c", "#7b3f9e",
    "#1f77b4", "#d62728", "#bcbd22", "#17becf", "#8c564b",
    "#e377c2", "#2ca02c", "#ff7f0e", "#9467bd", "#aec7e8",
    "#ffbb78", "#9e9e9e", "#cccccc",  # Localization, Other
]


def _compile_categories():
    out = []
    for name, spec in FUNCTIONAL_CATEGORIES:
        gos = set(spec.get("go_ids", ()))
        kws = [_re.compile(p, _re.IGNORECASE) for p in spec.get("keywords", ())]
        out.append((name, gos, kws))
    return out


def classify_proteins(flat_rows):
    """Map each protein to all categories its GO terms touch.

    Returns a Counter of `category_name -> protein_membership_count`.
    A protein in 3 functional categories contributes 1 to each (so the sum
    is the total number of category memberships, not the number of proteins).

    Fallback handling for proteins matching no functional category:
        * If at least one of their GO terms is a subcellular-localization
          term (cytoplasm, nucleus, mitochondrion, Golgi, ...), tag them
          "Subcellular localization only" -- we know where they live, just
          not what they do.
        * Otherwise tag them "Other / unclassified".
    This split keeps the residual `Other` slice genuinely uncharacterised.
    """
    by_protein = defaultdict(list)
    for acc, go_id, desc, ns, ri in flat_rows:
        by_protein[acc].append((go_id, desc))

    cats = _compile_categories()
    counter = Counter()
    for acc, terms in by_protein.items():
        matched = set()
        has_localization = False
        for go_id, desc in terms:
            if desc and _LOCALIZATION_KW.search(desc):
                has_localization = True
            for name, gos, kws in cats:
                if name in matched:
                    continue
                if go_id in gos:
                    matched.add(name)
                    continue
                if desc and any(rx.search(desc) for rx in kws):
                    matched.add(name)
        if matched:
            for name in matched:
                counter[name] += 1
        elif has_localization:
            counter[_LOC_ONLY_LABEL] += 1
        else:
            counter[_OTHER_LABEL] += 1
    return counter


def parse_results(path, min_score):
    """Stream results.csv: collect counters AND per-row tuples for the TSV."""
    proteins_total = set()
    proteins_high_conf = set()
    namespace_counts_all = Counter()
    namespace_counts_hc = Counter()
    total_rows = 0
    high_conf_rows = 0
    # Per-row tuples kept for the user-facing flat TSV. Only rows above the
    # cutoff are emitted (raw results.csv stays available alongside for users
    # who want the full unfiltered set).
    flat_rows = []

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_rows += 1
            acc = (row.get("query_accession") or row.get("accession") or "").strip()
            go_id = (row.get("go_id") or "").strip()
            ns = (row.get("category") or row.get("go_namespace") or "").strip()
            desc = (row.get("go_description") or row.get("go_name") or "").strip()
            try:
                ri = float(row.get("reliability_index") or row.get("max_ri") or 0)
            except ValueError:
                ri = 0.0
            if acc:
                proteins_total.add(acc)
            if ns:
                namespace_counts_all[ns] += 1
            if ri >= min_score:
                high_conf_rows += 1
                if acc:
                    proteins_high_conf.add(acc)
                if ns:
                    namespace_counts_hc[ns] += 1
                if acc and go_id:
                    flat_rows.append((acc, go_id, desc, ns, ri))

    return {
        "total_rows": total_rows,
        "high_conf_rows": high_conf_rows,
        "proteins_total": len(proteins_total),
        "proteins_high_conf": len(proteins_high_conf),
        "namespace_counts_all": namespace_counts_all,
        "namespace_counts_hc": namespace_counts_hc,
        "flat_rows": flat_rows,
    }


def write_go_terms_tsv(stats, out_path):
    """Flat per-(transcript, GO term) TSV with human-readable GO names."""
    rows = stats["flat_rows"]
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, "w") as f:
        f.write("transcript_id\tgo_id\tgo_name\tgo_namespace\treliability_index\n")
        for acc, go_id, desc, ns, ri in rows:
            f.write(f"{acc}\t{go_id}\t{desc}\t{ns}\t{ri:.4f}\n")


def write_summary(stats, min_score, out_path):
    with open(out_path, "w") as f:
        f.write("FANTASIA-Lite functional annotation summary\n")
        f.write("===========================================\n\n")
        f.write(f"Reliability-index cutoff: {min_score}\n\n")
        f.write(f"Total GO assignments:               {stats['total_rows']:,}\n")
        f.write(f"  ... above cutoff:                 {stats['high_conf_rows']:,}\n")
        f.write(f"Proteins with at least one GO term: {stats['proteins_total']:,}\n")
        f.write(f"  ... above cutoff:                 {stats['proteins_high_conf']:,}\n\n")
        f.write("GO assignments per namespace (above cutoff):\n")
        if not stats["namespace_counts_hc"]:
            f.write("  (none)\n")
        else:
            for ns, _ in sorted(
                stats["namespace_counts_hc"].items(),
                key=lambda kv: kv[1],
                reverse=True,
            ):
                label = NAMESPACE_LABELS.get(ns, ns)
                f.write(f"  {label}: {stats['namespace_counts_hc'][ns]:,}\n")


def write_categories_pie(stats, min_score, out_path):
    """Pie chart of broad functional categories.

    Each protein is counted once per matching curated category, so the pie
    represents fraction of category memberships rather than fraction of
    proteins. Proteins matching no category fall into "Other / unclassified".
    """
    counter = classify_proteins(stats["flat_rows"])
    if not counter:
        return False

    # Stable order: curated functional categories in their declared order,
    # then the "Localization only" bucket, then the "Other" bucket pinned
    # to the very end. Drop empty buckets.
    ordered = []
    for name, _ in FUNCTIONAL_CATEGORIES:
        if counter.get(name, 0) > 0:
            ordered.append((name, counter[name]))
    if counter.get(_LOC_ONLY_LABEL, 0) > 0:
        ordered.append((_LOC_ONLY_LABEL, counter[_LOC_ONLY_LABEL]))
    if counter.get(_OTHER_LABEL, 0) > 0:
        ordered.append((_OTHER_LABEL, counter[_OTHER_LABEL]))

    labels = [name for name, _ in ordered]
    values = [v for _, v in ordered]
    colors = _CATEGORY_COLORS[: len(labels)]
    total = sum(values)

    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    wedges, _texts, autotexts = ax.pie(
        values,
        labels=None,
        colors=colors,
        startangle=90,
        counterclock=False,
        # Only label slices >=2.5% of the pie to keep the figure readable
        autopct=lambda p: f"{p:.1f}%" if p >= 2.5 else "",
        pctdistance=0.78,
        wedgeprops={"edgecolor": "white", "linewidth": 1.0},
    )
    for at in autotexts:
        at.set_fontsize(10)
        at.set_color("white")
        at.set_fontweight("bold")
    legend_labels = [
        f"{name}  ({count:,}; {100 * count / total:.1f}%)"
        for name, count in ordered
    ]
    ax.legend(
        wedges,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=10,
        frameon=False,
        title=f"Category memberships  (RI \u2265 {min_score:.2f})",
        title_fontsize=11,
    )
    ax.set_title(
        "FANTASIA-Lite functional categories\n"
        "Each protein contributes once per matching curated category",
        fontsize=12,
        pad=14,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", required=True, help="Path to FANTASIA-Lite results.csv")
    p.add_argument("--out-dir", required=True, help="Directory to write summary + plot")
    p.add_argument("--min-score", type=float, default=0.5,
                   help="Reliability index cutoff for the high-confidence summary (default: 0.5)")
    args = p.parse_args()

    if not os.path.isfile(args.results):
        sys.stderr.write(f"results.csv not found: {args.results}\n")
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)
    stats = parse_results(args.results, args.min_score)

    summary_path = os.path.join(args.out_dir, "fantasia_summary.txt")
    plot_path = os.path.join(args.out_dir, "fantasia_go_categories.png")
    tsv_path = os.path.join(args.out_dir, "fantasia_go_terms.tsv")

    write_summary(stats, args.min_score, summary_path)
    write_categories_pie(stats, args.min_score, plot_path)
    write_go_terms_tsv(stats, tsv_path)

    print(f"Wrote {summary_path}")
    print(f"Wrote {plot_path}")
    print(f"Wrote {tsv_path}")


if __name__ == "__main__":
    main()
