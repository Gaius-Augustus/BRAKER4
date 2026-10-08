#!/usr/bin/env python3
"""
Write the lncRNA transcripts of FEELnc_codpot.pl as GFF3 lnc_RNA -> exon lines.

FEELnc writes <candidates>.lncRNA.gtf with exon lines only (no transcript
lines); the span of a transcript is the union of its exons. The transcripts
are numbered in genome order as <stem>-lncRNA_<n> (Name = ID, the StringTie
transcript ID as Alias, so that the transcript can be looked up in
feelnc_classifier.txt), exons as <id>.exon<k>. merge_ncrna_gff3.py adds the
gene line. Comment lines of the GTF (the notes of rule run_feelnc on why
no lncRNA was called) are copied into the GFF3 header.

Usage:
    python3 feelnc_to_gff3.py --stem <sample> lncRNAs.gtf -o lncRNAs.gff3
"""

import argparse
import re
import sys

TX_ID = re.compile(r'transcript_id "([^"]+)"')


def read_exons(path):
    """({transcript_id: [(seqid, source, start, end, strand), ...]}, comment lines) in file order."""
    exons, comments = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            m = TX_ID.search(cols[8])
            if not m:
                sys.exit(f"feelnc_to_gff3.py: exon without transcript_id: {line.rstrip()}")
            exons.setdefault(m.group(1), []).append((cols[0], cols[1], int(cols[3]), int(cols[4]), cols[6]))
    return exons, comments


def to_gff3(exons, stem):
    """GFF3 lines (without header), one lnc_RNA and its exons per transcript."""
    transcripts = []
    for tid, ex in exons.items():
        seqid, source, strand = ex[0][0], ex[0][1], ex[0][4]
        if any(e[0] != seqid or e[4] != strand for e in ex):
            sys.exit(f"feelnc_to_gff3.py: exons of {tid} on several sequences or strands")
        transcripts.append((seqid, min(e[2] for e in ex), max(e[3] for e in ex), strand, source, tid, ex))
    transcripts.sort(key=lambda t: (t[0], t[1], t[2], t[5]))
    lines = []
    for n, (seqid, start, end, strand, source, tid, ex) in enumerate(transcripts, 1):
        rid = f"{stem}-lncRNA_{n}"
        lines.append("\t".join([seqid, source, "lnc_RNA", str(start), str(end), ".", strand, ".",
                                f"ID={rid};Name={rid};Alias={tid};biotype=lncRNA"]))
        for k, (_, _, es, ee, _) in enumerate(sorted(ex, key=lambda e: e[2]), 1):
            lines.append("\t".join([seqid, source, "exon", str(es), str(ee), ".", strand, ".",
                                    f"ID={rid}.exon{k};Parent={rid}"]))
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gtf", help="<candidates>.lncRNA.gtf of FEELnc_codpot.pl (exon lines)")
    ap.add_argument("--stem", required=True, help="ID prefix: <stem>-lncRNA_<n>")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    exons, comments = read_exons(args.gtf)
    with open(args.output, "w") as out:
        out.write("##gff-version 3\n")
        for line in comments:
            out.write(line + "\n")
        for line in to_gff3(exons, args.stem):
            out.write(line + "\n")
    print(f"feelnc_to_gff3.py: {len(exons)} lncRNA transcripts", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
