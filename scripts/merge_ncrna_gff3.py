#!/usr/bin/env python3
"""
Merge ncRNA annotations into a protein-coding GFF3.

The protein-coding GFF3 is copied unchanged. The ncRNA GFF3 files follow in
priority order (e.g. rRNA, tRNA, Infernal, lncRNA). Every ncRNA is written as
gene -> RNA -> exon(s), with gene_biotype on the gene, as in NCBI/Ensembl GFF3
(needed for Annotrieve and other tools that count genes by biotype):

    X1  tRNAscan-SE  gene  121982  122066  .     -  .  ID=gene-s-X1.trna1;gene_biotype=tRNA
    X1  tRNAscan-SE  tRNA  121982  122066  64.2  -  .  ID=s-X1.trna1;Parent=gene-s-X1.trna1;...
    X1  tRNAscan-SE  exon  122029  122066  .     -  .  ID=s-X1.trna1.exon1;Parent=s-X1.trna1

An ncRNA gene is dropped when more than --max-coding-overlap of its exon
length overlaps protein-coding exons/CDS on the same strand, or when more
than --max-ncrna-overlap of its span overlaps an ncRNA gene on the same
strand from a file earlier in the priority list (e.g. an Infernal tRNA hit
at a locus that tRNAscan-SE already annotated).

Features that are not RNAs (e.g. Infernal hits to cis-regulatory Rfam
families such as riboswitches) are appended unchanged, without a gene.

Usage:
    merge_ncrna_gff3.py --coding braker.gff3 \
        --ncrna rRNA.gff3 tRNAs.gff3 ncRNAs_infernal.gff3 lncRNAs.gff3 \
        -o braker_with_ncRNA.gff3
"""

import argparse
import os
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
RNA_TYPES = {
    "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "scaRNA", "miRNA",
    "pre_miRNA", "lnc_RNA", "antisense_RNA", "ribozyme", "RNase_P_RNA",
    "RNase_MRP_RNA", "SRP_RNA", "telomerase_RNA", "vault_RNA", "Y_RNA",
    "transcript",
}
# gene_biotype for an RNA feature type when the RNA carries no biotype attribute
BIOTYPE_OF_TYPE = {"lnc_RNA": "lncRNA", "transcript": "ncRNA"}


def parse_attrs(col9):
    attrs = {}
    for part in col9.strip().strip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def format_attrs(attrs):
    first = [k for k in ("ID", "Parent") if k in attrs]
    rest = [k for k in attrs if k not in ("ID", "Parent")]
    return ";".join(f"{k}={attrs[k]}" for k in first + rest)


def read_gff3(path):
    """Return (header_lines, features); a feature is a list of 9 columns."""
    header, features = [], []
    if not path or not os.path.exists(path):
        return header, features
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                header.append(line)
                continue
            cols = line.split("\t")
            if len(cols) == 9:
                features.append(cols)
    return header, features


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


class IntervalIndex:
    """Intervals per (seqid, strand), kept sorted, disjoint and merged.

    Queries find the overlapping intervals by bisection, so one query costs
    O(log n + hits) instead of a scan over every stored interval (a scan made
    the merge quadratic: 50k coding intervals x 20k ncRNAs took two minutes).
    Added intervals are buffered and folded in on the next query: a large batch
    (the coding exons) is sorted and merged in one pass, a few intervals at a
    time (the accepted ncRNAs) are inserted by bisection.
    """

    BULK = 16  # pending intervals above which a flush re-sorts everything

    def __init__(self):
        self.starts = defaultdict(list)   # merged intervals, sorted by start
        self.ends = defaultdict(list)     # ends[k] belongs to starts[k]
        self.pending = defaultdict(list)  # added since the last query

    def add(self, seqid, strand, start, end):
        self.pending[(seqid, strand)].append((start, end))

    def _flush(self, key):
        pending = self.pending.pop(key, None)
        if not pending:
            return
        starts, ends = self.starts[key], self.ends[key]
        if len(pending) > self.BULK:
            merged_s, merged_e = [], []
            for s, e in sorted(pending + list(zip(starts, ends))):
                if merged_e and s <= merged_e[-1] + 1:
                    if e > merged_e[-1]:
                        merged_e[-1] = e
                else:
                    merged_s.append(s)
                    merged_e.append(e)
            self.starts[key], self.ends[key] = merged_s, merged_e
            return
        for s, e in pending:
            # intervals [i, j) overlap or touch [s, e]; replace them by the union
            i = bisect_left(ends, s - 1)
            j = bisect_right(starts, e + 1)
            if i < j:
                s, e = min(s, starts[i]), max(e, ends[j - 1])
            starts[i:j] = [s]
            ends[i:j] = [e]

    def covered_bp(self, seqid, strand, start, end):
        """bp of [start, end] covered by the stored intervals (merged)."""
        key = (seqid, strand)
        self._flush(key)
        starts, ends = self.starts.get(key), self.ends.get(key)
        if not starts:
            return 0
        i = bisect_left(ends, start)      # first interval ending at or after start
        j = bisect_right(starts, end)     # first interval starting after end
        return sum(min(ends[k], end) - max(starts[k], start) + 1 for k in range(i, j))


def build_loci(features, label, log):
    """Group features into loci: a root feature and all its descendants."""
    by_id, children, roots = {}, defaultdict(list), []
    for f in features:
        fid = parse_attrs(f[8]).get("ID")
        if fid:
            by_id[fid] = f
    for f in features:
        parents = parse_attrs(f[8]).get("Parent")
        if not parents:
            roots.append(f)
            continue
        parent = parents.split(",")[0]
        if parent in by_id:
            children[parent].append(f)
        else:
            log.append(f"[WARN] {label}: dropped {f[2]} {f[0]}:{f[3]}-{f[4]}, "
                       f"Parent={parent} not in file")
    loci = []
    for root in roots:
        members, stack = [root], [root]
        while stack:
            fid = parse_attrs(stack.pop()[8]).get("ID")
            for child in children.get(fid, ()):
                members.append(child)
                stack.append(child)
        loci.append(members)
    return loci


def rna_type_of(rna):
    """SO type of an RNA line; FEELnc lncRNAs come as transcript + biotype=lncRNA."""
    if rna[2] == "lncRNA" or (rna[2] == "transcript"
                              and parse_attrs(rna[8]).get("biotype") == "lncRNA"):
        return "lnc_RNA"
    return rna[2]


def biotype_of(rna):
    attrs = parse_attrs(rna[8])
    return (attrs.get("gene_biotype") or attrs.get("biotype")
            or BIOTYPE_OF_TYPE.get(rna[2], rna[2]))


def normalize_locus(members):
    """Return the locus as gene -> RNA -> exon lines, or None if not an RNA."""
    root = members[0]
    if root[2] in GENE_TYPES:
        rnas = [m[:2] + [rna_type_of(m)] + m[3:] for m in members[1:]
                if m[2] in RNA_TYPES or m[2] == "lncRNA"]
        if not rnas:
            return None
        gene_attrs = parse_attrs(root[8])
        gene_attrs.setdefault("gene_biotype", biotype_of(rnas[0]))
        gene = root[:8] + [format_attrs(gene_attrs)]
    elif root[2] in RNA_TYPES or root[2] == "lncRNA":
        rna_attrs = parse_attrs(root[8])
        rna = root[:2] + [rna_type_of(root)] + root[3:]
        gene_attrs = {"ID": f"gene-{rna_attrs['ID']}", "gene_biotype": biotype_of(rna)}
        if "Name" in rna_attrs:
            gene_attrs["Name"] = rna_attrs["Name"]
        rna_attrs["Parent"] = gene_attrs["ID"]
        rnas = [rna[:8] + [format_attrs(rna_attrs)]]
        gene = root[:2] + ["gene"] + root[3:5] + [".", root[6], ".", format_attrs(gene_attrs)]
    else:
        return None

    out = [gene]
    for rna in rnas:
        out.append(rna)
        rid = parse_attrs(rna[8])["ID"]
        exons = [m for m in members if m[2] == "exon"
                 and parse_attrs(m[8]).get("Parent", "").split(",")[0] == rid]
        if not exons:  # RNAs without exons get one exon spanning the RNA
            exons = [rna[:2] + ["exon"] + rna[3:5] + [".", rna[6], ".",
                                                     f"ID={rid}.exon1;Parent={rid}"]]
        out.extend(exons)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--coding", required=True, help="protein-coding GFF3 (copied unchanged)")
    ap.add_argument("--ncrna", nargs="*", default=[],
                    help="ncRNA GFF3 files, highest priority first; missing files are skipped")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--max-coding-overlap", type=float, default=0.5)
    ap.add_argument("--max-ncrna-overlap", type=float, default=0.5)
    args = ap.parse_args()

    log = []
    header, coding = read_gff3(args.coding)
    coding_idx = IntervalIndex()
    for f in coding:
        if f[2] in ("exon", "CDS"):
            coding_idx.add(f[0], f[6], int(f[3]), int(f[4]))

    accepted_idx = IntervalIndex()
    appended = []
    for path in args.ncrna:
        label = os.path.basename(path)
        if not os.path.exists(path):
            log.append(f"[INFO] {label}: not present, skipped")
            continue
        _, features = read_gff3(path)
        n_genes = n_other = n_coding = n_dup = 0
        for members in build_loci(features, label, log):
            locus = normalize_locus(members)
            if locus is None:
                appended.extend(members)
                n_other += 1
                continue
            gene = locus[0]
            seqid, strand, start, end = gene[0], gene[6], int(gene[3]), int(gene[4])
            exon_bp = sum(int(e[4]) - int(e[3]) + 1 for e in locus if e[2] == "exon")
            coding_bp = sum(coding_idx.covered_bp(seqid, strand, int(e[3]), int(e[4]))
                            for e in locus if e[2] == "exon")
            if exon_bp and coding_bp / exon_bp > args.max_coding_overlap:
                n_coding += 1
                continue
            if accepted_idx.covered_bp(seqid, strand, start, end) / (end - start + 1) \
                    > args.max_ncrna_overlap:
                n_dup += 1
                continue
            accepted_idx.add(seqid, strand, start, end)
            appended.extend(locus)
            n_genes += 1
        log.append(f"[INFO] {label}: {n_genes} ncRNA genes added, {n_other} other features "
                   f"added, {n_coding} dropped (overlap with coding exons), {n_dup} dropped "
                   f"(overlap with higher-priority ncRNA)")

    with open(args.output, "w") as out:
        if not any(h.startswith("##gff-version") for h in header):
            out.write("##gff-version 3\n")
        for h in header:
            out.write(h + "\n")
        for f in coding:
            out.write("\t".join(f) + "\n")
        for f in appended:
            out.write("\t".join(f) + "\n")
    print("\n".join(log), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
