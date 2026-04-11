#!/usr/bin/env python3
"""
Convert Infernal cmscan --tblout output to GFF3 format.

Parses the tabular output from cmscan (--fmt 2 --tblout) and produces
GFF3 features for ncRNA hits that pass Rfam gathering thresholds.

Usage:
    python3 infernal_to_gff3.py -i input.tblout -o output.gff3 [-p sample_prefix]
"""

import argparse
import sys


def parse_tblout(tblout_file):
    """Parse Infernal --fmt 2 tblout output.

    Fields (fmt 2):
    0: idx, 1: target name, 2: accession, 3: query name, 4: accession,
    5: clan name, 6: mdl, 7: mdl from, 8: mdl to, 9: seq from, 10: seq to,
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

            target_name = fields[0]  # sequence name
            rfam_acc = fields[3]     # e.g. RF00001
            rfam_name = fields[2]    # e.g. 5S_rRNA
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
                'seqid': target_name,
                'source': 'Infernal',
                'type': 'ncRNA',
                'start': start,
                'end': end,
                'score': score,
                'strand': gff_strand,
                'evalue': evalue,
                'rfam_acc': rfam_acc,
                'rfam_name': rfam_name,
            })

    return hits


def write_gff3(hits, output_file, prefix=''):
    """Write hits as GFF3."""
    with open(output_file, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, hit in enumerate(hits, 1):
            feature_id = f"{prefix}ncRNA_{i}" if prefix else f"ncRNA_{i}"
            attrs = (
                f"ID={feature_id};"
                f"Name={hit['rfam_name']};"
                f"Dbxref=RFAM:{hit['rfam_acc']};"
                f"evalue={hit['evalue']};"
                f"note=Infernal cmscan hit to {hit['rfam_name']} ({hit['rfam_acc']})"
            )
            fh.write(
                f"{hit['seqid']}\t{hit['source']}\t{hit['type']}\t"
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
    args = parser.parse_args()

    hits = parse_tblout(args.input)
    write_gff3(hits, args.output, prefix=args.prefix + '-' if args.prefix else '')

    print(f"Converted {len(hits)} Infernal hits to GFF3", file=sys.stderr)


if __name__ == '__main__':
    main()
