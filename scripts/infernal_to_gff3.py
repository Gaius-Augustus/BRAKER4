#!/usr/bin/env python3
"""
Convert Infernal cmscan --tblout output to GFF3 format.

Parses the tabular output from cmscan (--fmt 2 --tblout) and produces
GFF3 features for ncRNA hits that pass Rfam gathering thresholds.

With --family-types (rfam_family_types.tsv, from Rfam family.txt), each hit
gets the RNA type of its Rfam family (tRNA, snoRNA, miRNA, ...) and a
gene_biotype; merge_ncrna_gff3.py then adds the gene feature. Hits to
cis-regulatory families and self-splicing introns are not RNA genes and get
an SO type without a gene_biotype. Without the table every hit is an ncRNA.

Usage:
    python3 infernal_to_gff3.py -i input.tblout -o output.gff3 [-p sample_prefix] \
        [--family-types rfam_family_types.tsv]
"""

import argparse
import os
import sys


def load_family_types(path):
    """Rfam accession -> Rfam type string (e.g. 'Gene; snRNA; snoRNA; CD-box;')."""
    types = {}
    if not path or not os.path.exists(path):
        return types
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("rfam_acc\t"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                types[cols[0]] = cols[2]
    return types


def rna_class(rfam_type):
    """(GFF3 feature type, gene_biotype or None) for an Rfam type string.

    Feature types and biotypes follow NCBI/Ensembl GFF3. gene_biotype None
    means the hit is not an RNA gene (cis-regulatory element, intron).
    """
    parts = [p.strip() for p in rfam_type.split(";") if p.strip()]
    if not parts:
        return "ncRNA", "ncRNA"
    if parts[0] == "Cis-reg":
        if "riboswitch" in parts:
            return "riboswitch", None
        if "IRES" in parts:
            return "internal_ribosome_entry_site", None
        return "regulatory_region", None
    if parts[0] == "Intron":
        return "autocatalytically_spliced_intron", None
    for key, feature, biotype in (
        ("tRNA", "tRNA", "tRNA"),
        ("rRNA", "rRNA", "rRNA"),
        ("scaRNA", "scaRNA", "scaRNA"),
        ("snoRNA", "snoRNA", "snoRNA"),
        ("snRNA", "snRNA", "snRNA"),
        ("miRNA", "miRNA", "miRNA"),
        ("lncRNA", "lnc_RNA", "lncRNA"),
        ("antisense", "antisense_RNA", "antisense_RNA"),
        ("ribozyme", "ribozyme", "ribozyme"),
    ):
        if key in parts:
            return feature, biotype
    return "ncRNA", "ncRNA"


def parse_tblout(tblout_file):
    """Parse Infernal cmscan --fmt 2 tblout output.

    For cmscan, "target" is the Rfam CM model and "query" is the input
    genome sequence — so the GFF3 seqid must come from the query name.

    Fields (fmt 2):
    0: idx, 1: target name (Rfam family), 2: accession (RFxxxxx),
    3: query name (genome seq), 4: accession, 5: clan name,
    6: mdl, 7: mdl from, 8: mdl to, 9: seq from, 10: seq to,
    11: strand, 12: trunc, 13: pass, 14: gc, 15: bias, 16: score,
    17: E-value, 18: inc, 19: olp, 20: anyidx, 21: apts1, 22: apts2,
    23: winidx, 24: wpts1, 25: wpts2, 26: description of target
    """
    hits = []
    with open(tblout_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split()
            if len(fields) < 27:
                continue

            # Only keep included hits (inc == '!')
            if fields[18] != '!':
                continue

            # Only keep non-overlapping winners (olp == '*')
            if fields[19] != '*':
                continue

            seqid = fields[3]        # query name = genome sequence (scaffold)
            rfam_name = fields[1]    # e.g. 5S_rRNA
            rfam_acc = fields[2]     # e.g. RF00001
            seq_from = int(fields[9])
            seq_to = int(fields[10])
            strand = fields[11]
            score = fields[16]
            evalue = fields[17]

            # Normalize coordinates (GFF3 is 1-based, start < end)
            if seq_from > seq_to:
                start, end = seq_to, seq_from
            else:
                start, end = seq_from, seq_to

            gff_strand = '+' if strand == '+' else '-'

            hits.append({
                'seqid': seqid,
                'source': 'Infernal',
                'start': start,
                'end': end,
                'score': score,
                'strand': gff_strand,
                'evalue': evalue,
                'rfam_acc': rfam_acc,
                'rfam_name': rfam_name,
            })

    return hits


def write_gff3(hits, output_file, prefix='', family_types=None):
    """Write hits as GFF3."""
    family_types = family_types or {}
    with open(output_file, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, hit in enumerate(hits, 1):
            feature_id = f"{prefix}ncRNA_{i}" if prefix else f"ncRNA_{i}"
            rfam_type = family_types.get(hit['rfam_acc'], "")
            feature_type, biotype = rna_class(rfam_type)
            attrs = (
                f"ID={feature_id};"
                f"Name={hit['rfam_name']};"
                f"Dbxref=RFAM:{hit['rfam_acc']};"
                f"evalue={hit['evalue']};"
                f"note=Infernal cmscan hit to {hit['rfam_name']} ({hit['rfam_acc']})"
            )
            if biotype:
                attrs += f";gene_biotype={biotype}"
            if rfam_type:
                # ';' separates GFF3 attributes, so the type parts are comma-joined
                attrs += ";rfam_type=" + ",".join(
                    p.strip() for p in rfam_type.split(";") if p.strip())
            fh.write(
                f"{hit['seqid']}\t{hit['source']}\t{feature_type}\t"
                f"{hit['start']}\t{hit['end']}\t{hit['score']}\t"
                f"{hit['strand']}\t.\t{attrs}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description='Convert Infernal cmscan tblout (--fmt 2) to GFF3'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Infernal tblout file (--fmt 2)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output GFF3 file')
    parser.add_argument('-p', '--prefix', default='',
                        help='Prefix for feature IDs (e.g. sample name)')
    parser.add_argument('--family-types', default='',
                        help='rfam_family_types.tsv (Rfam accession -> type)')
    args = parser.parse_args()

    family_types = load_family_types(args.family_types)
    if args.family_types and not family_types:
        print(f"WARNING: no Rfam family types read from {args.family_types}; "
              "all hits are written as ncRNA", file=sys.stderr)
    hits = parse_tblout(args.input)
    write_gff3(hits, args.output, prefix=args.prefix + '-' if args.prefix else '',
               family_types=family_types)

    print(f"Converted {len(hits)} Infernal hits to GFF3", file=sys.stderr)


if __name__ == '__main__':
    main()
