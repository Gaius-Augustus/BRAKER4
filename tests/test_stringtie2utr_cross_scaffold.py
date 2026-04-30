"""
Regression tests for the cross-scaffold MSTRG ID collision bug.

Root cause (issue #29): When StringTie is run per-scaffold and outputs are
merged, MSTRG transcript IDs restart at 1 per scaffold and collide. read_gtf()
groups all features under the same transcript_id key, so a single entry in
stringtie_gtf can contain exons from multiple scaffolds. If a BRAKER multi-exon
transcript on scaffold_1 is matched to that entry (correctly, via shared
scaffold_1 introns), merge_features() was adding scaffold_3 exons too, causing
fix_feature_coordinates() to compute wildly wrong gene bounds.

Fix: merge_features() now filters StringTie exons to the same seqname as the
BRAKER transcript before any further processing. fix_feature_coordinates() also
skips features whose seqname does not match the transcript line as a safety net.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stringtie2utr import merge_features, fix_feature_coordinates  # noqa: E402


def _gtf_line(seqname, source, feature, start, end, score, strand, frame, attrs):
    return "\t".join([seqname, source, feature,
                      str(start), str(end),
                      score, strand, frame, attrs])


# ---------------------------------------------------------------------------
# Test 1 – REGRESSION: multi-exon BRAKER on scaffold_1 matched to a StringTie
#   entry that has exons on BOTH scaffold_1 (correct) and scaffold_3 (collision).
#   Expected: scaffold_3 exons must NOT appear; gene bounds stay within scaffold_1.
# ---------------------------------------------------------------------------

def test_cross_scaffold_mstrg_collision_multi_exon():
    """Exons from a colliding scaffold must not be merged into the BRAKER transcript."""

    braker_tx_id = "g1545.t2"
    stringtie_tx_id = "MSTRG.1.1"

    # BRAKER multi-exon transcript on scaffold_1 (1000-3000)
    braker_features = [
        _gtf_line("scaffold_1", "AUGUSTUS", "exon", 1000, 1500, ".", "+", ".",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
        _gtf_line("scaffold_1", "AUGUSTUS", "CDS",  1000, 1497, ".", "+", "0",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
        _gtf_line("scaffold_1", "AUGUSTUS", "exon", 2000, 3000, ".", "+", ".",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
        _gtf_line("scaffold_1", "AUGUSTUS", "CDS",  2000, 2997, ".", "+", "0",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
    ]

    # StringTie entry for MSTRG.1.1 — contains exons from two scaffolds because
    # per-scaffold StringTie runs produced colliding MSTRG.1.1 IDs.
    stringtie_features = [
        # correct scaffold_1 exons (match the BRAKER intron pattern)
        _gtf_line("scaffold_1", "StringTie", "exon",  900, 1500, "1000", "+", ".",
                  'gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "1";'),
        _gtf_line("scaffold_1", "StringTie", "exon", 2000, 3200, "1000", "+", ".",
                  'gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "2";'),
        # wrong scaffold_3 exons from a colliding MSTRG.1.1 on another scaffold
        _gtf_line("scaffold_3", "StringTie", "exon", 6763, 8867, "1000", "+", ".",
                  'gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "1";'),
    ]

    tsebra_gtf = {braker_tx_id: list(braker_features)}
    stringtie_gtf = {stringtie_tx_id: list(stringtie_features)}
    selected = {braker_tx_id: stringtie_tx_id}

    result = merge_features(tsebra_gtf, stringtie_gtf, selected)
    merged = result[braker_tx_id]

    seqnames = [line.split('\t')[0] for line in merged]
    starts   = [int(line.split('\t')[3]) for line in merged]
    ends     = [int(line.split('\t')[4]) for line in merged]

    assert "scaffold_3" not in seqnames, (
        "scaffold_3 features must not be merged into a scaffold_1 BRAKER transcript"
    )
    assert 8867 not in ends, "scaffold_3 exon end (8867) must not appear"
    assert 6763 not in starts, "scaffold_3 exon start (6763) must not appear"
    assert max(ends) == 3200, f"Expected max end 3200 (scaffold_1 UTR extension), got {max(ends)}"
    assert min(starts) == 900, f"Expected min start 900 (scaffold_1 UTR extension), got {min(starts)}"


# ---------------------------------------------------------------------------
# Test 2 – REGRESSION: fix_feature_coordinates must not use cross-scaffold
#   features when computing transcript/gene coordinate bounds.
# ---------------------------------------------------------------------------

def test_fix_feature_coordinates_ignores_wrong_scaffold():
    """fix_feature_coordinates must skip features whose seqname differs from the transcript."""

    tx_id = "g1545.t2"
    gene_id_attr = 'gene_id "g1545";'

    # Transcript line on scaffold_1, originally spanning 1000-3000
    tx_line = _gtf_line("scaffold_1", "AUGUSTUS", "transcript",
                         1000, 3000, ".", "+", ".", 'transcript_id "g1545.t2"; gene_id "g1545";')
    gene_line = _gtf_line("scaffold_1", "AUGUSTUS", "gene",
                           1000, 3000, ".", "+", ".", gene_id_attr)

    tx_dict   = {tx_id: tx_line}
    gene_dict = {gene_id_attr: gene_line}
    tx_to_gene_dict = {tx_id: gene_id_attr}

    # Features: two scaffold_1 exons + one rogue scaffold_3 feature
    gtf_dict = {tx_id: [
        _gtf_line("scaffold_1", "AUGUSTUS", "exon", 1000, 1500, ".", "+", ".",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
        _gtf_line("scaffold_1", "AUGUSTUS", "exon", 2000, 3200, ".", "+", ".",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
        # wrong-scaffold feature that must be ignored
        _gtf_line("scaffold_3", "stringtie2utr", "five_prime_UTR", 6763, 8867, ".", "+", ".",
                  'transcript_id "g1545.t2"; gene_id "g1545";'),
    ]}

    gene_dict_out, tx_dict_out = fix_feature_coordinates(
        gtf_dict, gene_dict, tx_to_gene_dict, tx_dict)

    tx_fields   = tx_dict_out[tx_id].split('\t')
    gene_fields = gene_dict_out[gene_id_attr].split('\t')

    assert int(tx_fields[3]) == 1000, f"Transcript start must be 1000, got {tx_fields[3]}"
    assert int(tx_fields[4]) == 3200, f"Transcript end must be 3200, got {tx_fields[4]}"
    assert int(gene_fields[3]) == 1000, f"Gene start must be 1000, got {gene_fields[3]}"
    assert int(gene_fields[4]) == 3200, f"Gene end must be 3200, got {gene_fields[4]}"
