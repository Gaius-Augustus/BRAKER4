#!/usr/bin/env python3
"""
Write a BRAKER4 GFF3 as GFF3 for ENA decoupled annotation (Webin-CLI, -context
annotation, ANALYSIS_TYPE DECOUPLED_ANNOTATION).

Format: https://ena-docs.readthedocs.io/en/latest/submit/gff3.html
Rules followed here are the ones ENA's validator gff3tools
(https://github.com/enasequence/gff3tools) enforces, as of its main branch
of 2026-11-16:

  * ##gff-version 3, optional ##species with the NCBI taxonomy URI, and one
    ##sequence-region before the features of each sequence.
  * Every attribute except ID, Parent and partial becomes an INSDC qualifier,
    so only INSDC qualifiers are written: locus_tag, product, transl_table,
    codon_start, inference, partial, pseudo, note, Dbxref.
  * Feature types must map to an INSDC feature: start_codon, stop_codon and
    intron rows are dropped; RNA types are kept (they map to ncRNA with an
    ncRNA_class, or to tRNA, rRNA).
  * A feature on several rows (spliced CDS, UTR, RNA) repeats one ID; rows are
    in ascending order and do not overlap; Parent comes before the children.
  * Phase is set on CDS rows only and recomputed from the CDS. The reading
    frame of a CDS that does not start at phase 0 is given as codon_start,
    which is what the INSDC translation reads.
  * locus_tag: PREFIX_number, PREFIX 3-12 upper-case letters/digits starting
    with a letter (registered with ENA), on the gene and all its features.
  * partial: a model without start codon (or that does not start with one) is
    5' partial, one without in-frame stop codon is 3' partial. gff3tools reads
    partial in genomic orientation (start = lower coordinate, end = upper
    coordinate, also on the minus strand: GFF3Feature.isFivePrimePartial and
    the complement test fixtures), so that is what is written. The partial
    attribute goes on the terminal row(s) of the CDS, and on mRNA, exon and
    gene where they end at the same coordinate.
  * Complete CDS shorter than 25 aa are rejected without evidence, so every
    CDS carries inference="ab initio prediction:<program>[:<version>]".
  * A CDS with an internal stop codon is written as pseudo (what ENA's own
    translation fix would do) and reported.
  * tRNA: product tRNA-Xxx; complete tRNAs must be 50-150 bp (summed exons),
    must not lie inside the CDS of a gene; rRNA must not overlap a CDS or
    another rRNA on the same strand. Features breaking these rules are left
    out and reported.

The genome FASTA gives the sequence lengths for ##sequence-region and the
codons for the start/stop/frame checks. Sequences are read one at a time, so
memory stays at about one chromosome.

Usage:
    gff3_to_ena.py --gff3 braker.gff3 --genome genome.fa \\
        --locus-tag-prefix BRAKR --taxid 4932 --translation-table 1 \\
        [--seqid-map seqids.tsv] [--software-versions software_versions.tsv] \\
        -o braker.ena.gff3 [--id-map braker.ena.ids.tsv] [--report report.txt]
"""

import argparse
import re
import sys
from collections import Counter, defaultdict

# ---------------------------------------------------------------------------
# Genetic codes as gff3tools defines them (TranslationTableDescriptor):
# amino acids and start codons in TCAG order. BRAKER4 supports 1, 6 and 29.
# ---------------------------------------------------------------------------
GENETIC_CODES = {
    1: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M---------------M---------------M----------------------------"),
    6: ("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------"),
    29: ("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "--------------*--------------------M----------------------------"),
}

LOCUS_TAG_PREFIX_RE = re.compile(r"^[A-Z][A-Z0-9]{2,11}$")
# gff3tools parses ##sequence-region with (?<accessionId>[^.]+)(\.\d+)?
ENA_SEQID_RE = re.compile(r"^[^.\s]+(\.\d+)?$")
SPECIES_URI = "https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{}"

GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
TRANSCRIPT_TYPES = {"mRNA", "transcript"}
# RNA types written as they are; gff3tools maps each to an INSDC feature
# (feature-mapping.tsv): tRNA, rRNA, or ncRNA with /ncRNA_class.
RNA_TYPES = {
    "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "scaRNA", "miRNA", "pre_miRNA",
    "lnc_RNA", "antisense_RNA", "ribozyme", "RNase_P_RNA", "RNase_MRP_RNA",
    "SRP_RNA", "telomerase_RNA", "vault_RNA", "Y_RNA",
    "autocatalytically_spliced_intron",
}
RNA_ALIASES = {"lncRNA": "lnc_RNA", "pseudogenic_tRNA": "tRNA",
               "pseudogenic_rRNA": "rRNA"}
# Features without gene that are written on their own (Infernal hits to
# cis-regulatory Rfam families).
STANDALONE_TYPES = {"riboswitch": "riboswitch",
                    "regulatory_region": "regulatory_region",
                    "internal_ribosome_entry_site": "regulatory_region"}
# Mature miRNAs are about 22 nt; Rfam miRNA families model the hairpin.
MATURE_MIRNA_MAX_LEN = 30
TRNA_MIN_LEN, TRNA_MAX_LEN = 50, 150

AMINO_ACIDS = {
    "Ala": "Ala", "Arg": "Arg", "Asn": "Asn", "Asp": "Asp", "Cys": "Cys",
    "Gln": "Gln", "Glu": "Glu", "Gly": "Gly", "His": "His", "Ile": "Ile",
    "Leu": "Leu", "Lys": "Lys", "Met": "Met", "Phe": "Phe", "Pro": "Pro",
    "Ser": "Ser", "Thr": "Thr", "Trp": "Trp", "Tyr": "Tyr", "Val": "Val",
    "iMet": "Met", "fMet": "fMet", "SeC": "Sec", "Sec": "Sec", "Pyl": "Pyl",
}

RESERVED = {"%": "%25", ";": "%3B", "=": "%3D", "&": "%26", ",": "%2C",
            # gff3tools decodes values with URLDecoder, which turns '+' into
            # a space.
            "+": "%2B"}


class ConversionError(Exception):
    pass


def encode(value):
    value = re.sub(r"[\t\r\n]+", " ", str(value)).strip()
    value = "".join(RESERVED.get(ch, ch) for ch in value)
    if not value.isascii():
        value = value.encode("ascii", "replace").decode("ascii")
    return value


def decode(value):
    return re.sub(r"%([0-9A-Fa-f]{2})", lambda m: chr(int(m.group(1), 16)), value)


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

class Feature:
    __slots__ = ("seqid", "source", "type", "start", "end", "score", "strand",
                 "phase", "attrs", "line")

    def __init__(self, cols, line):
        self.seqid, self.source, self.type = cols[0], cols[1], cols[2]
        self.start, self.end = int(cols[3]), int(cols[4])
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        self.score, self.strand, self.phase = cols[5], cols[6], cols[7]
        self.attrs = {}
        for part in cols[8].strip().strip(";").split(";"):
            if "=" in part:
                key, value = part.split("=", 1)
                self.attrs[key.strip()] = decode(value.strip())
        self.line = line

    @property
    def id(self):
        return self.attrs.get("ID")

    @property
    def parents(self):
        p = self.attrs.get("Parent")
        return [x for x in p.split(",") if x] if p else []


def read_gff3(path):
    features = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise ConversionError(f"{path}:{n}: expected 9 columns, got {len(cols)}")
            features.append(Feature(cols, n))
    return features


def read_fasta(path):
    """Yield (seqid, sequence) one record at a time."""
    seqid, chunks = None, []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if seqid is not None:
                    yield seqid, "".join(chunks).upper()
                seqid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if seqid is not None:
        yield seqid, "".join(chunks).upper()


def read_seqid_map(path):
    mapping = {}
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                raise ConversionError(f"{path}:{n}: expected <seqid>TAB<ENA seqid>")
            mapping[cols[0]] = cols[1]
    return mapping


def read_inference(path):
    """Program name -> version string from software_versions.tsv."""
    programs = {}
    if not path:
        return programs
    try:
        fh = open(path)
    except OSError:
        return programs
    with fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            name, version = cols[0].strip(), cols[1].strip()
            m = re.match(r"v?(\d+(?:\.\d+)+)\b", version)
            programs.setdefault(name, m.group(1) if m else "")
    return programs


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

COMPLEMENT = str.maketrans("ACGTRYKMBVDHNSW", "TGCAYRMKVBHDNSW")
BASES = "TCAG"


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


class Code:
    def __init__(self, table):
        if table not in GENETIC_CODES:
            raise ConversionError(f"translation table {table} not supported "
                                  f"(supported: {sorted(GENETIC_CODES)})")
        aas, starts = GENETIC_CODES[table]
        self.table = table
        self.other, self.start = {}, {}
        i = 0
        for b1 in BASES:
            for b2 in BASES:
                for b3 in BASES:
                    codon = b1 + b2 + b3
                    self.other[codon] = aas[i]
                    self.start[codon] = "M" if starts[i] == "M" else aas[i]
                    i += 1

    def is_start(self, codon):
        return self.start.get(codon) == "M"

    def is_stop(self, codon):
        return self.other.get(codon) == "*"

    def translate(self, seq):
        return "".join(self.other.get(seq[i:i + 3], "X")
                       for i in range(0, len(seq) - 2, 3))


def merge_intervals(intervals):
    """Sorted, merged (start, end) list; abutting intervals are joined."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def subtract(intervals, remove):
    """Parts of `intervals` not covered by `remove` (both lists of (start, end))."""
    out = []
    for s, e in merge_intervals(intervals):
        pieces = [(s, e)]
        for rs, re_ in remove:
            nxt = []
            for ps, pe in pieces:
                if re_ < ps or rs > pe:
                    nxt.append((ps, pe))
                    continue
                if ps < rs:
                    nxt.append((ps, rs - 1))
                if pe > re_:
                    nxt.append((re_ + 1, pe))
            pieces = nxt
        out.extend(pieces)
    return out


def phases(segments_bio, phase0):
    """Phase of each CDS segment, segments in biological order."""
    result, p = [], phase0
    for s, e in segments_bio:
        result.append(p)
        rest = (e - s + 1 - p) % 3
        p = (3 - rest) % 3
    return result


# ---------------------------------------------------------------------------
# Model building
# ---------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.counts = Counter()
        self.dropped_types = Counter()
        self.messages = []

    def note(self, msg):
        self.messages.append(msg)


class OutFeature:
    """One output feature; rows are (start, end) segments sharing an ID."""

    def __init__(self, ftype, source, strand, rows, fid=None, parent=None):
        self.type, self.source, self.strand = ftype, source, strand
        self.rows = sorted(rows)
        self.id, self.parent = fid, parent
        self.attrs = {}          # INSDC qualifiers, written in insertion order
        self.row_phases = None   # CDS only, in ascending row order
        self.partial_lo = False  # partial at the lower coordinate
        self.partial_hi = False  # partial at the upper coordinate

    @property
    def start(self):
        return self.rows[0][0]

    @property
    def end(self):
        return self.rows[-1][1]


class Gene:
    def __init__(self, seqid, gid, source, strand, kind):
        self.seqid, self.id, self.source, self.strand = seqid, gid, source, strand
        self.kind = kind            # "coding" or "rna" or "standalone"
        self.gene = None            # OutFeature
        self.children = []          # list of lists: transcript + its features
        self.braker_ids = []        # (gene_id, transcript_id) for the id map

    @property
    def start(self):
        return min(f.start for group in self.children for f in group)

    @property
    def end(self):
        return max(f.end for group in self.children for f in group)


def program_inference(source, programs):
    """INSDC /inference for a gene model predicted by `source`."""
    if "augustus" in source.lower():
        name = "AUGUSTUS"
    elif source.lower() == "gmst":
        # GeneMarkS-T, run by GeneMark-ETP on assembled transcripts
        name = "GeneMarkS-T"
    elif "genemark" in source.lower():
        name = next((p for p in programs if p.startswith("GeneMark")), "GeneMark")
    else:
        name = source
    version = programs.get(name, "")
    name = name.replace(" ", "_").replace(":", "_")
    return f"ab initio prediction:{name}" + (f":{version}" if version else "")


def build_coding_transcript(tx, kids, seq, code, programs, report):
    """Return the output features of one protein-coding transcript, or None."""
    strand = tx.strand
    if strand not in "+-":
        report.note(f"{tx.id}: strand '{strand}' is not + or -, transcript left out")
        report.counts["dropped_transcript_no_strand"] += 1
        return None
    cds_in = [k for k in kids if k.type == "CDS"]
    stops = [k for k in kids if k.type == "stop_codon"]
    starts = [k for k in kids if k.type == "start_codon"]
    exons_in = [(k.start, k.end) for k in kids if k.type == "exon"]
    # CDS includes the stop codon in INSDC; merge stop_codon rows in case the
    # source wrote it outside the CDS.
    cds = merge_intervals([(k.start, k.end) for k in cds_in] +
                          [(k.start, k.end) for k in stops])
    # Without start_codon and stop_codon rows (GeneMark transcripts in some
    # braker.gtf files) completeness comes from the sequence, and a stop codon
    # right after the CDS, inside the transcript's exon, belongs to the CDS.
    infer = not starts and not stops
    bio = cds if strand == "+" else cds[::-1]
    # The reading frame comes from the phase of the 5'-most CDS row; the stop
    # codon is only ever merged at the 3' end, so it does not shift the frame.
    first_cds = min(cds_in, key=lambda k: k.start) if strand == "+" \
        else max(cds_in, key=lambda k: k.end)
    phase0 = int(first_cds.phase) if first_cds.phase in ("0", "1", "2") else 0

    nt ="".join(seq[s - 1:e] for s, e in cds)
    if strand == "-":
        nt = revcomp(nt)
    if len(nt) != sum(e - s + 1 for s, e in cds):
        raise ConversionError(f"{tx.id}: CDS beyond the end of sequence {tx.seqid}")

    if infer and (len(nt) - phase0) % 3 == 0 and not code.is_stop(nt[-3:]):
        span = exons_in or [(tx.start, tx.end)]
        nxt = (cds[-1][1] + 1, cds[-1][1] + 3) if strand == "+" else \
            (cds[0][0] - 3, cds[0][0] - 1)
        inside = nxt[0] >= 1 and any(s <= nxt[0] and nxt[1] <= e for s, e in span)
        codon = seq[nxt[0] - 1:nxt[1]] if inside else ""
        if inside and code.is_stop(codon if strand == "+" else revcomp(codon)):
            cds = merge_intervals(cds + [nxt])
            bio = cds if strand == "+" else cds[::-1]
            nt = nt + (codon if strand == "+" else revcomp(codon))
            report.counts["stop_codon_taken_from_sequence"] += 1

    five_complete = ((bool(starts) or infer) and phase0 == 0 and len(nt) >= 3
                     and code.is_start(nt[:3]))
    if starts and not five_complete:
        report.note(f"{tx.id}: has start_codon but CDS does not begin with a "
                    f"start codon of table {code.table} in frame; written 5' partial")
    coding = nt[phase0:]
    whole = len(coding) - len(coding) % 3
    three_complete = ((bool(stops) or infer) and len(coding) % 3 == 0 and whole >= 3
                      and code.is_stop(coding[whole - 3:whole]))
    if infer:
        report.counts["completeness_from_sequence"] += 1
    if stops and not three_complete:
        report.note(f"{tx.id}: has stop_codon but CDS does not end with an "
                    f"in-frame stop codon of table {code.table}; written 3' partial")
    protein = code.translate(coding[:whole])
    body = protein[:-1] if three_complete else protein
    pseudo = "*" in body
    if pseudo:
        report.note(f"{tx.id}: internal stop codon at amino acid "
                    f"{body.index('*') + 1}; CDS written as pseudo")
        report.counts["pseudo_internal_stop"] += 1
    if five_complete and three_complete and not pseudo and len(body) < 25:
        report.counts["complete_cds_lt_25aa"] += 1

    # exons; UTRs never overlap the CDS (stringtie2utr can label a stop codon
    # that is missing from the CDS as 3' UTR)
    utr5_in = subtract([(k.start, k.end) for k in kids if k.type == "five_prime_UTR"], cds)
    utr3_in = subtract([(k.start, k.end) for k in kids if k.type == "three_prime_UTR"], cds)
    exons = merge_intervals(exons_in) if exons_in else []
    if not exons or not all(any(es <= s and e <= ee for es, ee in exons) for s, e in cds):
        if exons:
            report.note(f"{tx.id}: exons do not contain the CDS; exons rebuilt from CDS and UTRs")
        exons = merge_intervals(cds + utr5_in + utr3_in + exons_in)
    # A UTR lies beyond a start or stop codon. On the side where the CDS is
    # partial there is none: the transcript ends with the CDS (AGAT turns an
    # AUGUSTUS intron that runs off the sequence end into a UTR).
    for complete, utr_rows, side_lo in ((five_complete, utr5_in, strand == "+"),
                                        (three_complete, utr3_in, strand == "-")):
        if complete:
            continue
        if utr_rows:
            report.counts["utr_removed_at_partial_end"] += 1
            utr_rows.clear()
        if side_lo:
            exons = [(max(s, cds[0][0]), e) for s, e in exons if e >= cds[0][0]]
        else:
            exons = [(s, min(e, cds[-1][1])) for s, e in exons if s <= cds[-1][1]]
    lo, hi = exons[0][0], exons[-1][1]

    src = tx.source
    feats = []
    mrna = OutFeature("mRNA", src, strand, [(lo, hi)], fid=tx.id)
    feats.append(mrna)
    for s, e in exons:
        feats.append(OutFeature("exon", src, strand, [(s, e)], parent=tx.id))
    cds_f = OutFeature("CDS", src, strand, cds, fid=f"cds-{tx.id}", parent=tx.id)
    bio_phases = phases(bio, phase0)
    cds_f.row_phases = bio_phases if strand == "+" else bio_phases[::-1]
    cds_f.attrs["product"] = "hypothetical protein"
    cds_f.attrs["inference"] = program_inference(src, programs)
    if pseudo:
        cds_f.attrs.pop("product")
        cds_f.attrs["pseudo"] = "true"
        cds_f.attrs["note"] = "internal stop codon"
    else:
        cds_f.attrs["transl_table"] = str(code.table)
    if phase0:
        cds_f.attrs["codon_start"] = str(phase0 + 1)
    feats.append(cds_f)
    for ftype, rows, tag in (("five_prime_UTR", utr5_in, "utr5"),
                             ("three_prime_UTR", utr3_in, "utr3")):
        rows = merge_intervals(rows)
        if rows:
            feats.append(OutFeature(ftype, src, strand, rows,
                                    fid=f"{tag}-{tx.id}", parent=tx.id))

    # partial: 5' end is the lower coordinate on +, the upper on -
    five_lo = strand == "+"
    for complete, at_lo in ((five_complete, five_lo), (three_complete, not five_lo)):
        if complete:
            continue
        report.counts["partial_5prime" if at_lo == five_lo else "partial_3prime"] += 1
        coord = cds[0][0] if at_lo else cds[-1][1]
        for f in feats:
            edge = f.start if at_lo else f.end
            if f.type == "CDS" or (f.type in ("mRNA", "exon") and edge == coord):
                if at_lo:
                    f.partial_lo = True
                else:
                    f.partial_hi = True
    return feats


def rna_product_and_note(rna, ftype, report):
    attrs = rna.attrs
    product, notes = None, []
    if ftype == "tRNA":
        iso = attrs.get("isotype")
        if not iso:
            m = re.search(r"tRNA\d*-([A-Za-z]+?)([ACGTU]{3})?$", attrs.get("Name", ""))
            iso = m.group(1) if m else ""
        aa = AMINO_ACIDS.get(iso)
        product = f"tRNA-{aa}" if aa else "tRNA-OTHER"
        if iso == "iMet":
            notes.append("initiator tRNA")
        if attrs.get("anticodon"):
            notes.append(f"anticodon {attrs['anticodon']}")
    elif ftype == "rRNA":
        product = attrs.get("product", "")
        product = re.sub(r"\s*\(partial\)\s*$", "", product)
        if not product:
            name = attrs.get("Name", "")
            for pat, prod in ((r"5[._]8S", "5.8S ribosomal RNA"),
                              (r"(^|[^.\d])5S", "5S ribosomal RNA"),
                              (r"18S|SSU", "18S ribosomal RNA"),
                              (r"25S", "25S ribosomal RNA"),
                              (r"28S|LSU", "28S ribosomal RNA")):
                if re.search(pat, name):
                    product = prod
                    break
        product = product or None
    for key in ("note", "Note"):
        if attrs.get(key):
            notes.append(attrs[key])
    dbx = [x for x in attrs.get("Dbxref", "").split(",") if x]
    return product, notes, dbx


def build_rna_transcript(tx, kids, report):
    ftype = RNA_ALIASES.get(tx.type, tx.type)
    if ftype == "transcript" and tx.attrs.get("biotype") == "lncRNA":
        ftype = "lnc_RNA"
    if ftype not in RNA_TYPES:
        return None, None
    rows = merge_intervals([(k.start, k.end) for k in kids if k.type == "exon"]) \
        or [(tx.start, tx.end)]
    if ftype == "miRNA" and sum(e - s + 1 for s, e in rows) > MATURE_MIRNA_MAX_LEN:
        ftype = "pre_miRNA"
    pseudo = tx.type.startswith("pseudogenic_") or tx.attrs.get("pseudo") == "true"
    rna = OutFeature(ftype, tx.source, tx.strand, rows, fid=tx.id)
    product, notes, dbx = rna_product_and_note(tx, ftype, report)
    if pseudo:
        rna.attrs["pseudo"] = "true"
        if product:
            notes.insert(0, f"{product} pseudogene")
    elif product:
        rna.attrs["product"] = product
    if notes:
        rna.attrs["note"] = "; ".join(notes)
    if dbx:
        rna.attrs["Dbxref"] = dbx
    return rna, pseudo


def index_features(features):
    """ID index, children by parent ID, and features without parent."""
    by_id = {}
    children = defaultdict(list)
    for f in features:
        if f.id:
            if f.id in by_id and (f.type in GENE_TYPES or f.type in TRANSCRIPT_TYPES):
                raise ConversionError(f"line {f.line}: ID {f.id} used twice")
            by_id.setdefault(f.id, f)
    for f in features:
        for p in f.parents:
            children[p].append(f)
    roots = [f for f in features if not f.parents]
    return by_id, children, roots


# ---------------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------------

def convert(args):
    report = Report()
    code = Code(args.translation_table)
    prefix = args.locus_tag_prefix.strip()
    if not LOCUS_TAG_PREFIX_RE.match(prefix):
        raise ConversionError(
            f"locus_tag prefix '{prefix}' is invalid: ENA needs 3-12 upper-case "
            f"letters or digits, starting with a letter (e.g. BRAKR)")
    programs = read_inference(args.software_versions)
    seqid_map = read_seqid_map(args.seqid_map) if args.seqid_map else {}

    features = read_gff3(args.gff3)
    by_id, children, roots = index_features(features)

    # transcripts of each gene; orphan transcripts get a gene
    genes_by_seq = defaultdict(list)
    used_ids = set(by_id)

    def unique_id(base):
        fid, n = base, 1
        while fid in used_ids:
            n += 1
            fid = f"{base}-{n}"
        used_ids.add(fid)
        return fid

    root_items = []
    for r in roots:
        if r.type in GENE_TYPES:
            txs = [c for c in children.get(r.id, []) if c.type not in GENE_TYPES]
            root_items.append((r, txs))
        elif r.type in TRANSCRIPT_TYPES or r.type in RNA_TYPES or r.type in RNA_ALIASES:
            root_items.append((None, [r]))
        elif r.type in STANDALONE_TYPES:
            root_items.append((r, None))
        else:
            report.dropped_types[r.type] += 1

    seqids_with_features = []
    for gene_f, txs in root_items:
        anchor = gene_f or txs[0]
        seqid = anchor.seqid
        if seqid not in genes_by_seq:
            seqids_with_features.append(seqid)
        if txs is None:   # standalone regulatory feature
            g = Gene(seqid, anchor.id, anchor.source, anchor.strand, "standalone")
            ftype = STANDALONE_TYPES[anchor.type]
            out = OutFeature(ftype, anchor.source, anchor.strand,
                             [(anchor.start, anchor.end)], fid=None)
            _, notes, dbx = rna_product_and_note(anchor, ftype, report)
            if anchor.type == "internal_ribosome_entry_site":
                notes.insert(0, "internal ribosome entry site")
            if anchor.type == "riboswitch" and not notes:
                notes.append("riboswitch")
            if notes:
                out.attrs["note"] = "; ".join(notes)
            if dbx:
                out.attrs["Dbxref"] = dbx
            if not out.attrs:
                out.attrs["note"] = anchor.type.replace("_", " ")
            g.children.append([out])
            genes_by_seq[seqid].append(g)
            report.counts["standalone_" + ftype] += 1
            continue
        g = Gene(seqid, gene_f.id if gene_f else None, anchor.source, anchor.strand, None)
        gene_pseudo = bool(gene_f) and (gene_f.type == "pseudogene"
                                        or gene_f.attrs.get("pseudo") == "true")
        pseudo_flags = []
        for tx in sorted(txs, key=lambda t: (t.start, t.end, t.id or "")):
            if tx.seqid != seqid:
                raise ConversionError(f"{tx.id}: on {tx.seqid}, its gene on {seqid}")
            kids = children.get(tx.id, [])
            has_cds = any(k.type == "CDS" for k in kids)
            if tx.type in TRANSCRIPT_TYPES and has_cds:
                if g.kind == "rna":
                    report.note(f"{tx.id}: coding transcript in an RNA gene, left out")
                    continue
                g.kind = "coding"
                g.children.append(("coding", tx, kids))
            else:
                rna, pseudo = build_rna_transcript(tx, kids, report)
                if rna is None:
                    report.dropped_types[f"{tx.type} (no CDS)"] += 1
                    continue
                if g.kind == "coding":
                    report.note(f"{tx.id}: non-coding transcript in a protein-coding "
                                f"gene, left out")
                    continue
                g.kind = "rna"
                g.children.append([rna])
                pseudo_flags.append(pseudo or gene_pseudo)
            g.braker_ids.append((gene_f.id if gene_f else "", tx.id))
        if not g.children:
            continue
        if g.id is None:
            g.id = unique_id(f"gene-{txs[0].id}")
        g.pseudo_rna = g.kind == "rna" and all(pseudo_flags) and bool(pseudo_flags)
        if g.kind == "rna" and gene_pseudo:
            for group in g.children:
                group[0].attrs["pseudo"] = "true"
                group[0].attrs.pop("product", None)
        genes_by_seq[seqid].append(g)

    # ENA sequence names, checked before the genome is read
    ena_names = {}
    for seqid in seqids_with_features:
        name = seqid_map.get(seqid, seqid) if seqid_map else seqid
        if seqid_map and seqid not in seqid_map:
            raise ConversionError(f"sequence {seqid} has features but is not in "
                                  f"the seqid map {args.seqid_map}")
        if not ENA_SEQID_RE.match(name):
            raise ConversionError(
                f"sequence name '{name}' cannot be read by ENA (only letters, digits "
                f"and symbols without '.', optionally followed by .<version>); "
                f"give a seqid map to the ENA accessions")
        ena_names[seqid] = name
    if len(set(ena_names.values())) != len(ena_names):
        raise ConversionError("two sequences map to the same ENA name")

    # sequences: stream the genome, convert each sequence's genes
    out_blocks = {}
    seq_lengths = {}
    needed = set(genes_by_seq)
    for seqid, seq in read_fasta(args.genome):
        if seqid not in needed:
            continue
        seq_lengths[seqid] = len(seq)
        out_blocks[seqid] = finish_sequence(seqid, seq, genes_by_seq[seqid],
                                            code, programs, report)
        needed.discard(seqid)
    if needed:
        raise ConversionError("sequences with features missing from the genome: "
                              + ", ".join(sorted(needed)[:10]))

    # IDs must identify one feature in the whole file
    seen = set()
    for s in out_blocks:
        for g in out_blocks[s]:
            for f in ([g.gene] if g.gene else []) + [f for grp in g.children for f in grp]:
                if f.id:
                    if f.id in seen:
                        raise ConversionError(f"ID {f.id} would be used by two features")
                    seen.add(f.id)

    # locus tags in output order
    order = [s for s in seqids_with_features if out_blocks.get(s)]
    n_genes = sum(1 for s in order for g in out_blocks[s] if g.kind != "standalone")
    width = max(6, len(str(n_genes * 10)))
    counter = 0
    id_map = []
    for s in order:
        for g in out_blocks[s]:
            if g.kind == "standalone":
                continue
            counter += 10
            tag = f"{prefix}_{counter:0{width}d}"
            for group in [[g.gene]] + g.children:
                for f in group:
                    f.attrs = {"locus_tag": tag, **f.attrs}
            for gid, tid in g.braker_ids:
                id_map.append((ena_names[s], gid, tid, tag))

    with open(args.output, "w") as out:
        out.write("##gff-version 3\n")
        if args.taxid:
            out.write(f"##species {SPECIES_URI.format(args.taxid)}\n")
        for s in order:
            name = ena_names[s]
            out.write(f"##sequence-region {name} 1 {seq_lengths[s]}\n")
            for g in out_blocks[s]:
                for f in ([g.gene] if g.gene else []) + [f for grp in g.children for f in grp]:
                    write_feature(out, name, f)
    if args.id_map:
        with open(args.id_map, "w") as fh:
            fh.write("ena_seqid\tbraker_gene_id\tbraker_transcript_id\tlocus_tag\n")
            for row in id_map:
                fh.write("\t".join(row) + "\n")
    report.counts["genes"] = n_genes
    report.counts["sequences"] = len(order)
    return report


def finish_sequence(seqid, seq, genes, code, programs, report):
    """Build output features of all genes on one sequence; apply overlap rules."""
    built = []
    for g in genes:
        if g.kind == "coding":
            groups = []
            for _, tx, kids in g.children:
                feats = build_coding_transcript(tx, kids, seq, code, programs, report)
                if feats:
                    groups.append(feats)
                    report.counts["mrna"] += 1
            if not groups:
                continue
            g.children = groups
        for group in g.children:
            for f in group:
                if f.end > len(seq):
                    raise ConversionError(f"{f.id or f.parent}: ends at {f.end}, "
                                          f"beyond {seqid} ({len(seq)} bp)")
        if g.kind != "standalone":
            lo, hi = g.start, g.end
            gene = OutFeature("gene", g.source, g.strand, [(lo, hi)], fid=g.id)
            for group in g.children:
                t = group[0]
                t.parent = g.id
                gene.partial_lo |= t.partial_lo and t.start == lo
                gene.partial_hi |= t.partial_hi and t.end == hi
            if getattr(g, "pseudo_rna", False):
                gene.attrs["pseudo"] = "true"
            g.gene = gene
        built.append(g)

    # overlap rules for tRNA and rRNA against CDS (gff3tools FeatureOverlapValidation)
    cds_feats = [f for g in built if g.kind == "coding"
                 for group in g.children for f in group if f.type == "CDS"]
    cds_rows = [(f.strand, s, e) for f in cds_feats for s, e in f.rows]
    kept, rrnas = [], []
    for g in built:
        if g.kind == "rna":
            rna = g.children[0][0]
            if rna.type == "tRNA" and "pseudo" not in rna.attrs:
                length = sum(e - s + 1 for s, e in rna.rows)
                if not TRNA_MIN_LEN <= length <= TRNA_MAX_LEN:
                    report.note(f"{rna.id}: tRNA of {length} bp, ENA accepts "
                                f"{TRNA_MIN_LEN}-{TRNA_MAX_LEN} bp; left out")
                    report.counts["dropped_trna_length"] += 1
                    continue
            # every tRNA row inside a row of one CDS, either strand
            if rna.type == "tRNA" and any(
                    all(any(cs <= s and e <= ce for cs, ce in c.rows) for s, e in rna.rows)
                    for c in cds_feats if c.start <= rna.end and rna.start <= c.end):
                report.note(f"{rna.id}: tRNA inside a CDS; left out")
                report.counts["dropped_trna_in_cds"] += 1
                continue
            if rna.type == "rRNA":
                hit = any(st == rna.strand and s <= ce and cs <= e
                          for st, cs, ce in cds_rows for s, e in rna.rows)
                if hit:
                    report.note(f"{rna.id}: rRNA overlaps a CDS on the same strand; left out")
                    report.counts["dropped_rrna_overlap_cds"] += 1
                    continue
                if any(o.strand == rna.strand and o.start <= rna.end and rna.start <= o.end
                       for o in rrnas):
                    report.note(f"{rna.id}: rRNA overlaps another rRNA; left out")
                    report.counts["dropped_rrna_overlap_rrna"] += 1
                    continue
                rrnas.append(rna)
            report.counts["rna_" + rna.type] += 1
        kept.append(g)
    kept.sort(key=lambda g: (g.start, g.end, g.id or ""))
    for g in kept:
        g.children.sort(key=lambda grp: (grp[0].start, grp[0].end, grp[0].id or ""))
        for grp in g.children:
            head, rest = grp[0], grp[1:]
            rest.sort(key=lambda f: (f.start, f.end, f.type))
            grp[:] = [head] + rest
    return kept


def write_feature(out, seqid, f):
    """Write one feature, one row per segment; values are percent-encoded,
    list values (several Dbxref, partial=start,end) are joined by a literal
    comma, which gff3tools reads as a value separator."""
    rows = f.rows
    base = []
    if f.id:
        base.append(f"ID={encode(f.id)}")
    if f.parent:
        base.append(f"Parent={encode(f.parent)}")
    for key, value in f.attrs.items():
        values = value if isinstance(value, list) else [value]
        base.append(f"{key}=" + ",".join(encode(v) for v in values))
    for i, (s, e) in enumerate(rows):
        attrs = list(base)
        partial = []
        if f.partial_lo and i == 0:
            partial.append("start")
        if f.partial_hi and i == len(rows) - 1:
            partial.append("end")
        if partial:
            attrs.append("partial=" + ",".join(partial))
        col9 = ";".join(attrs)
        phase = str(f.row_phases[i]) if f.type == "CDS" else "."
        out.write("\t".join([seqid, f.source or ".", f.type, str(s), str(e), ".",
                             f.strand if f.strand in "+-" else ".", phase, col9]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--gff3", required=True, help="BRAKER4 GFF3 (braker.gff3 or braker_with_ncRNA.gff3)")
    ap.add_argument("--genome", required=True, help="genome FASTA the annotation refers to")
    ap.add_argument("--locus-tag-prefix", required=True,
                    help="locus_tag prefix registered with ENA")
    ap.add_argument("--taxid", default="", help="NCBI taxonomy ID for ##species")
    ap.add_argument("--translation-table", type=int, default=1)
    ap.add_argument("--seqid-map", default="",
                    help="TSV: BRAKER sequence name, ENA sequence name/accession")
    ap.add_argument("--software-versions", default="",
                    help="software_versions.tsv, for program versions in /inference")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--id-map", default="", help="TSV: BRAKER gene/transcript ID -> locus_tag")
    ap.add_argument("--report", default="", help="conversion report (default: stderr)")
    args = ap.parse_args()
    if args.taxid and not args.taxid.isdigit():
        sys.exit(f"ERROR: --taxid must be a number, got '{args.taxid}'")
    try:
        report = convert(args)
    except ConversionError as e:
        sys.exit(f"ERROR: {e}")
    lines = [f"{k}\t{v}" for k, v in sorted(report.counts.items())]
    lines += [f"dropped feature type\t{t}\t{n}" for t, n in sorted(report.dropped_types.items())]
    lines += report.messages
    text = "\n".join(lines) + "\n"
    if args.report:
        with open(args.report, "w") as fh:
            fh.write(text)
    else:
        sys.stderr.write(text)


if __name__ == "__main__":
    main()
