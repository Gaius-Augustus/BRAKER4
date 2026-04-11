#!/usr/bin/env python3
"""Classify missed and partially-correct reference transcripts in a BRAKER run
and check whether augustus.hints.gtf, genemark.gtf, or hintsfile.gff would
allow them to be rescued. Intended for the athaliana benchmark scenarios.

Usage:
    low_hanging_fruit_audit.py \
        --reference ref_annot.gff3 \
        --braker    braker.gtf \
        --augustus  augustus.hints.gtf \
        --genemark  genemark.gtf \
        --hints     hintsfile.gff \
        --label     ET
"""
import argparse
import re
import sys
from collections import defaultdict, Counter

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
ID_RE = re.compile(r'(?:ID|transcript_id)[ =]"?([^";]+)"?')
PARENT_RE = re.compile(r'Parent=([^;]+)')


def parse_attrs(field):
    d = {}
    for k, v in ATTR_RE.findall(field):
        d[k] = v
    if "=" in field and not d:
        for kv in field.split(";"):
            kv = kv.strip()
            if not kv or "=" not in kv:
                continue
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def load_cds_transcripts(path):
    """Return dict transcript_id -> dict(chrom, strand, cds=[(s,e)...]).
    Stop codons are folded into the adjacent CDS exon so that AUGUSTUS-style
    GTFs (stop_codon separate from CDS) compare equal to Phytozome-style
    references (stop included in CDS)."""
    tx = defaultdict(lambda: {"chrom": None, "strand": None,
                              "cds": [], "stops": []})
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype not in ("CDS", "stop_codon"):
                continue
            chrom, _, _, s, e, _, strand, _, attrs = f
            s = int(s); e = int(e)
            d = parse_attrs(attrs)
            tid = d.get("transcript_id")
            if tid is None:
                m = PARENT_RE.search(attrs)
                if m:
                    tid = m.group(1).split(",")[0]
            if tid is None:
                continue
            t = tx[tid]
            t["chrom"] = chrom
            t["strand"] = strand
            if ftype == "CDS":
                t["cds"].append((s, e))
            else:
                t["stops"].append((s, e))
    for t in tx.values():
        t["cds"].sort()
        for ss, se in t["stops"]:
            merged = False
            for i, (cs, ce) in enumerate(t["cds"]):
                if ce + 1 == ss:
                    t["cds"][i] = (cs, se)
                    merged = True
                    break
                if se + 1 == cs:
                    t["cds"][i] = (ss, ce)
                    merged = True
                    break
                if ss >= cs and se <= ce:
                    merged = True
                    break
            if not merged:
                t["cds"].append((ss, se))
        t["cds"].sort()
    return tx


def cds_summary(t):
    """Return (chrom, strand, intron_chain_tuple, cds_start, cds_end)."""
    cds = t["cds"]
    if not cds:
        return None
    introns = []
    for i in range(len(cds) - 1):
        introns.append((cds[i][1] + 1, cds[i + 1][0] - 1))
    return (t["chrom"], t["strand"], tuple(introns), cds[0][0], cds[-1][1])


def build_locus_index(transcripts):
    """Index transcripts on (chrom, strand) for fast overlap lookup."""
    idx = defaultdict(list)
    summaries = {}
    for tid, t in transcripts.items():
        s = cds_summary(t)
        if s is None:
            continue
        chrom, strand, introns, lo, hi = s
        idx[(chrom, strand)].append((lo, hi, tid, introns))
        summaries[tid] = s
    for k in idx:
        idx[k].sort()
    return idx, summaries


def overlapping(idx, chrom, strand, lo, hi):
    out = []
    for tlo, thi, tid, introns in idx.get((chrom, strand), []):
        if thi < lo:
            continue
        if tlo > hi:
            break
        out.append((tlo, thi, tid, introns))
    return out


def class_code(ref_summary, qry_summary):
    """Classify the prediction (qry) against a reference (ref) CDS transcript.

    Returns one of:
      '='                  : CDS-exact match (UTRs ignored)
      'bounds_only'        : same intron chain, only termini differ
      'extra_introns'      : qry has all ref introns plus extra ones (extra exons)
      'missing_introns'    : qry has fewer introns than ref (missing exons)
      'shared_junction'    : >=1 intron shared, neither subset
      'overlap_no_junction': same-strand overlap, no shared intron
      'no_overlap'         : nothing
    """
    rchrom, rstrand, rintrons, rlo, rhi = ref_summary
    qchrom, qstrand, qintrons, qlo, qhi = qry_summary
    if rchrom != qchrom or rstrand != qstrand:
        return "no_overlap"
    if qhi < rlo or qlo > rhi:
        return "no_overlap"
    if rintrons == qintrons and rlo == qlo and rhi == qhi:
        return "="
    rset = set(rintrons)
    qset = set(qintrons)
    shared = rset & qset
    if rintrons and qintrons and rset == qset:
        return "bounds_only"
    if rset and qset and rset.issubset(qset):
        return "extra_introns"
    if qset and rset and qset.issubset(rset):
        return "missing_introns"
    if shared:
        return "shared_junction"
    return "overlap_no_junction"


def best_match(ref_sum, hits, summaries):
    rank = {"=": 8, "bounds_only": 7, "missing_introns": 6, "extra_introns": 5,
            "shared_junction": 4, "overlap_no_junction": 3, "no_overlap": 1}
    best = ("no_overlap", None)
    for _, _, tid, _ in hits:
        c = class_code(ref_sum, summaries[tid])
        if rank[c] > rank[best[0]]:
            best = (c, tid)
    return best


def load_hint_intervals(path):
    """Return dict (chrom, strand, type) -> sorted list of (start, end, mult)."""
    h = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, ftype, s, e, _, strand, _, attrs = f
            if ftype not in ("intron", "start", "stop"):
                continue
            mult = 1
            m = re.search(r"mult=(\d+)", attrs)
            if m:
                mult = int(m.group(1))
            h[(chrom, strand, ftype)].append((int(s), int(e), mult))
    for k in h:
        h[k].sort()
    return h


def hint_lookup_set(hints, chrom, strand, ftype):
    return {(s, e) for s, e, _ in hints.get((chrom, strand, ftype), [])}


def hint_intron_coverage(ref_sum, intron_set):
    """Fraction of ref introns confirmed by hints set, exact-match."""
    _, _, introns, _, _ = ref_sum
    if not introns:
        return None
    matched = sum(1 for i in introns if i in intron_set)
    return matched / len(introns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--braker", required=True)
    ap.add_argument("--augustus", required=True)
    ap.add_argument("--genemark", required=True)
    ap.add_argument("--hints", required=True)
    ap.add_argument("--label", required=True)
    args = ap.parse_args()

    sys.stderr.write(f"[{args.label}] loading reference\n")
    ref_tx = load_cds_transcripts(args.reference)
    ref_idx, ref_sum = build_locus_index(ref_tx)
    sys.stderr.write(f"[{args.label}] loading braker\n")
    bra_tx = load_cds_transcripts(args.braker)
    bra_idx, bra_sum = build_locus_index(bra_tx)
    sys.stderr.write(f"[{args.label}] loading augustus\n")
    aug_tx = load_cds_transcripts(args.augustus)
    aug_idx, aug_sum = build_locus_index(aug_tx)
    sys.stderr.write(f"[{args.label}] loading genemark\n")
    gm_tx = load_cds_transcripts(args.genemark)
    gm_idx, gm_sum = build_locus_index(gm_tx)
    sys.stderr.write(f"[{args.label}] loading hints\n")
    hints = load_hint_intervals(args.hints)

    sys.stderr.write(f"[{args.label}] reference transcripts: {len(ref_sum)}\n")
    sys.stderr.write(f"[{args.label}] braker transcripts:    {len(bra_sum)}\n")
    sys.stderr.write(f"[{args.label}] augustus transcripts:  {len(aug_sum)}\n")
    sys.stderr.write(f"[{args.label}] genemark transcripts:  {len(gm_sum)}\n")

    bra_class_counter = Counter()
    rescue = Counter()
    rescue_pairs = Counter()
    intron_support_counter = Counter()
    sj_support = []  # list of (ref_class, augustus_class, genemark_class, intron_support)

    for tid, rsum in ref_sum.items():
        chrom, strand, introns, lo, hi = rsum
        bra_hits = overlapping(bra_idx, chrom, strand, lo, hi)
        bra_class, _ = best_match(rsum, bra_hits, bra_sum)
        bra_class_counter[bra_class] += 1
        if bra_class == "=":
            continue
        # not exact: try to rescue
        aug_hits = overlapping(aug_idx, chrom, strand, lo, hi)
        gm_hits = overlapping(gm_idx, chrom, strand, lo, hi)
        aug_class, _ = best_match(rsum, aug_hits, aug_sum)
        gm_class, _ = best_match(rsum, gm_hits, gm_sum)

        intron_set = hint_lookup_set(hints, chrom, strand, "intron")
        starts = hint_lookup_set(hints, chrom, strand, "start")
        stops = hint_lookup_set(hints, chrom, strand, "stop")
        cov = hint_intron_coverage(rsum, intron_set)
        # start/stop hint support of the ref CDS:
        rstart = (lo, lo + 2) if strand == "+" else (hi - 2, hi)
        rstop = (hi - 2, hi) if strand == "+" else (lo, lo + 2)
        has_start = rstart in starts
        has_stop = rstop in stops

        rescue[(bra_class, "aug=", aug_class == "=")] += 1
        rescue[(bra_class, "gm=", gm_class == "=")] += 1
        rescue_pairs[(bra_class, aug_class, gm_class)] += 1
        intron_support_counter[(bra_class,
                                "introns",
                                "all" if cov == 1.0 else
                                "some" if cov and cov > 0 else
                                "none" if introns else "single_exon",
                                "start_yes" if has_start else "start_no",
                                "stop_yes" if has_stop else "stop_no")] += 1

    print(f"\n========== {args.label}: braker class distribution (per ref transcript) ==========")
    total = sum(bra_class_counter.values())
    for c, n in bra_class_counter.most_common():
        print(f"  {c:24s} {n:7d}  {n/total*100:5.1f}%")
    print(f"  TOTAL                    {total:7d}")

    print(f"\n========== {args.label}: rescue eligibility for non-exact ref transcripts ==========")
    print("(rows = how braker classifies the ref; columns = does augustus / genemark have an EXACT CDS match for the ref?)")
    print(f"{'braker_class':>22s}  {'aug_exact':>10s}  {'gm_exact':>10s}  {'count':>8s}")
    pair_summary = Counter()
    for (bc, ac, gc), n in rescue_pairs.items():
        pair_summary[(bc, ac == "=", gc == "=")] += n
    for (bc, aug_eq, gm_eq), n in sorted(pair_summary.items(), key=lambda x: (-x[1])):
        print(f"{bc:>22s}  {str(aug_eq):>10s}  {str(gm_eq):>10s}  {n:8d}")

    print(f"\n========== {args.label}: detailed rescue pairs (top 30) ==========")
    print(f"{'braker_class':>22s}  {'aug_class':>22s}  {'gm_class':>22s}  {'count':>8s}")
    for (bc, ac, gc), n in rescue_pairs.most_common(30):
        print(f"{bc:>22s}  {ac:>22s}  {gc:>22s}  {n:8d}")

    print(f"\n========== {args.label}: hint support for non-exact ref transcripts ==========")
    print("(rows = braker_class | intron_support | start_hint | stop_hint)")
    by = Counter()
    for (bc, _, isup, sst, sto), n in intron_support_counter.items():
        by[(bc, isup, sst, sto)] += n
    for (bc, isup, sst, sto), n in sorted(by.items(), key=lambda x: (-x[1]))[:30]:
        print(f"  {bc:>22s} | introns={isup:>11s} | {sst:>9s} | {sto:>8s}  {n:7d}")


if __name__ == "__main__":
    main()
