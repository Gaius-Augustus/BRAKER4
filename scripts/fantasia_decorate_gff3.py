#!/usr/bin/env python3
"""
Decorate a BRAKER GFF3 file with GO term annotations from FANTASIA-Lite.

The Sequence Ontology / GFF3 specification reserves the `Ontology_term`
attribute for cross-references to controlled vocabularies. Multiple terms
are comma-separated, e.g.

    Ontology_term=GO:0008150,GO:0003674,GO:0005575

(see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
section "Reserved Attributes").

FANTASIA assigns GO terms at the protein level. Predicted-protein IDs match
the BRAKER transcript IDs that are encoded in the GFF3 -- but BRAKER4's GFF3
is unusual: it contains BOTH `transcript` features (AUGUSTUS / GeneMark
native, with `ID=g1.t1`) AND `mRNA` features (AGAT-converted, with
`ID=agat-mrna-1` and the original BRAKER ID stashed in `transcript_id=g1.t2`
as a separate attribute). To annotate both downstream conventions, this
script:

  * Decorates BOTH `mRNA` and `transcript` rows.
  * Uses `transcript_id=...` as the FANTASIA lookup key when present,
    falling back to `ID=...`. This makes AGAT-converted mRNA rows and
    AUGUSTUS/GeneMark-native transcript rows both resolve correctly.
  * Aggregates the per-transcript GO set to the parent gene feature so the
    gene also carries the rolled-up union.

Only assignments with reliability_index >= --min-score are kept. The original
GFF3 is preserved untouched; the decorated copy is written to --gff3-out.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


_ID_RE = re.compile(r'(?:^|;)\s*ID=([^;]+)')
_PARENT_RE = re.compile(r'(?:^|;)\s*Parent=([^;]+)')
_TXID_RE = re.compile(r'(?:^|;)\s*transcript_id=([^;]+)')
_ONT_RE = re.compile(r'(^|;)\s*Ontology_term=([^;]+)')

# BRAKER4 GFF3 has both `transcript` (native) and `mRNA` (AGAT) rows for
# essentially the same gene products. Decorate both so downstream tools
# that walk either feature type pick up the GO annotations.
_TRANSCRIPT_FEATURES = {"mRNA", "transcript"}


def _lookup_key(attr):
    """Return the FANTASIA lookup key for a transcript-equivalent row.

    Prefers `transcript_id=...` (set by AGAT on converted mRNA rows) and
    falls back to `ID=...` (used by AUGUSTUS/GeneMark-native transcript rows).
    """
    tx_m = _TXID_RE.search(attr)
    if tx_m:
        return tx_m.group(1)
    id_m = _ID_RE.search(attr)
    if id_m:
        return id_m.group(1)
    return None


def load_go_assignments(results_csv, min_score):
    """Map transcript_id -> set of GO IDs above the score threshold."""
    tx2go = defaultdict(set)
    with open(results_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tx = (row.get("query_accession") or row.get("accession") or "").strip()
            go = (row.get("go_id") or "").strip()
            try:
                ri = float(row.get("reliability_index") or row.get("max_ri") or 0)
            except ValueError:
                ri = 0.0
            if tx and go and ri >= min_score:
                tx2go[tx].add(go)
    return tx2go


def _add_ontology_term(attr, go_set):
    """Append or merge an Ontology_term attribute on a GFF3 col9 string.

    Preserves any existing Ontology_term value (defensively merging) so the
    decorator is idempotent and won't clobber annotations from other tools.
    """
    existing = _ONT_RE.search(attr)
    if existing:
        old_terms = {t for t in existing.group(2).split(",") if t}
        merged = sorted(old_terms.union(go_set))
        merged_value = ",".join(merged)
        # Replace the matched ";Ontology_term=..." (or "Ontology_term=..." at start)
        # in place so attribute order is preserved.
        leading = existing.group(1)
        return (
            attr[: existing.start()]
            + leading
            + f"Ontology_term={merged_value}"
            + attr[existing.end():]
        )
    new_value = ",".join(sorted(go_set))
    if not attr:
        return f"Ontology_term={new_value}"
    sep = "" if attr.endswith(";") else ";"
    return attr + sep + f"Ontology_term={new_value}"


def decorate(gff3_in, gff3_out, tx2go):
    """Two-pass: build transcript->gene map, then rewrite the GFF3.

    Both `mRNA` and `transcript` rows are decorated, and the lookup key for
    each row is `transcript_id=...` (preferred) or `ID=...` (fallback). See
    the module docstring for why BRAKER4's GFF3 needs this.
    """
    tx2gene = {}
    n_candidate_rows = 0
    with open(gff3_in) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in _TRANSCRIPT_FEATURES:
                continue
            n_candidate_rows += 1
            tx = _lookup_key(cols[8])
            par_m = _PARENT_RE.search(cols[8])
            if tx and par_m:
                tx2gene[tx] = par_m.group(1)

    gene2go = defaultdict(set)
    for tx, gos in tx2go.items():
        gene = tx2gene.get(tx)
        if gene:
            gene2go[gene].update(gos)

    n_tx = 0
    n_gene = 0
    with open(gff3_in) as fin, open(gff3_out, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(line)
                continue
            feat = cols[2]
            if feat in _TRANSCRIPT_FEATURES:
                tx = _lookup_key(cols[8])
                if tx and tx in tx2go:
                    cols[8] = _add_ontology_term(cols[8], tx2go[tx])
                    n_tx += 1
            elif feat == "gene":
                id_m = _ID_RE.search(cols[8])
                if id_m and id_m.group(1) in gene2go:
                    cols[8] = _add_ontology_term(cols[8], gene2go[id_m.group(1)])
                    n_gene += 1
            fout.write("\t".join(cols) + "\n")
    return n_tx, n_gene, n_candidate_rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gff3-in", required=True, help="Input GFF3 (e.g. braker.gff3)")
    p.add_argument("--gff3-out", required=True, help="Output decorated GFF3")
    p.add_argument("--results", required=True, help="FANTASIA-Lite results.csv")
    p.add_argument("--min-score", type=float, default=0.5,
                   help="Reliability index cutoff (default: 0.5)")
    args = p.parse_args()

    tx2go = load_go_assignments(args.results, args.min_score)
    n_tx, n_gene, n_candidates = decorate(args.gff3_in, args.gff3_out, tx2go)
    sys.stderr.write(
        f"Decorated {n_tx} transcript/mRNA and {n_gene} gene features with "
        f"Ontology_term GO IDs (cutoff RI >= {args.min_score}; "
        f"{n_candidates} candidate transcript-equivalent rows in input GFF3; "
        f"{len(tx2go)} transcript IDs in FANTASIA results)\n"
    )
    # Defensive: a fantasia run that produced GO assignments but couldn't
    # match a single transcript in the GFF3 is almost always an ID-format
    # mismatch -- fail loudly so the next pipeline run catches it instead of
    # silently emitting an unannotated GFF3.
    if n_candidates > 0 and len(tx2go) > 0 and n_tx == 0:
        sys.stderr.write(
            "ERROR: zero transcript features were decorated despite the input "
            "GFF3 containing transcript-equivalent rows AND the FANTASIA "
            "results.csv containing GO assignments. This is almost certainly "
            "an ID-format mismatch between braker.aa headers and the GFF3's "
            "transcript_id / ID attributes. Inspect a sample mRNA row and a "
            "sample query_accession from results.csv to debug.\n"
        )
        sys.exit(2)


if __name__ == "__main__":
    main()
