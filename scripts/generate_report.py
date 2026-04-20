#!/usr/bin/env python3
"""
Generate a self-contained HTML report for a BRAKER4 pipeline run.

The methods narrative is composed by `generate_methods_text(workdir, mode)`,
which inspects which output files exist on disk and writes one cohesive
paragraph in execution order. Rules do NOT contribute methods text directly.

Reads:
  - report_citations.txt: Plain-text citations (one per line, appended by rules)
  - report_citations.bib: BibTeX entries (appended by rules)
  - QC files: training_summary.png, busco_summary.txt, compleasm_summary.txt, etc.
  - Pipeline metadata: sample name, mode, date

Produces:
  - braker_report.html: Self-contained HTML with inline CSS and embedded images
  - braker_citations.bib: Deduplicated BibTeX file

Usage:
    python3 generate_report.py -d output/sample/ -o output/sample/results/
"""

import argparse
import base64
import os
import re
import sys
from collections import OrderedDict
from datetime import datetime



def read_file(path):
    """Read file contents, return empty string if not found."""
    if os.path.exists(path):
        with open(path) as f:
            return f.read()
    return ""


def embed_image(path, download_name=None):
    """Base64-encode an image for inline HTML embedding.

    If download_name is provided, a small download button is added below
    the image so users can save it directly from the HTML report.
    """
    if not os.path.exists(path):
        return ""
    with open(path, "rb") as f:
        raw = f.read()
    # Skip files that aren't actual images (e.g. error messages written as .png)
    if not raw.startswith(b'\x89PNG') and not raw.startswith(b'\xff\xd8') and not raw.startswith(b'<svg'):
        return ""
    data = base64.b64encode(raw).decode("ascii")
    ext = path.rsplit(".", 1)[-1].lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
            "svg": "image/svg+xml", "pdf": "application/pdf"}.get(ext, "image/png")
    data_uri = f"data:{mime};base64,{data}"
    img_tag = f'<img src="{data_uri}" style="max-width:100%; height:auto;">'
    if download_name:
        dl_btn = (
            f'<div style="text-align:right; margin-top:4px;">'
            f'<a download="{download_name}" href="{data_uri}" '
            f'style="font-size:0.8em; color:#3a7d44; text-decoration:none;">'
            f'\u2B73 Download {download_name}</a></div>'
        )
        return img_tag + dl_btn
    return img_tag


def read_software_versions(outdir, workdir):
    """Read software_versions.tsv and return list of (tool, version) tuples."""
    versions = OrderedDict()
    for search_dir in [outdir, workdir]:
        path = os.path.join(search_dir, "software_versions.tsv")
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t", 1)
                    if len(parts) == 2:
                        tool, version = parts
                        if tool not in versions:
                            versions[tool] = version
            break
    return list(versions.items())


def parse_bbc_log(path):
    """Parse a best_by_compleasm aggregated log into per-pass decisions.

    Returns a list of dicts (one per pass), each with keys:
        label, braker_missing, augustus_missing, genemark_missing,
        decision, new_missing
    Pass labels for dual mode are 'short-read', 'IsoSeq', 'merge';
    for non-dual modes there is a single unlabeled entry.
    Returns [] if the log is missing, empty, or marked as skipped.
    """
    if not os.path.exists(path):
        return []
    text = read_file(path).strip()
    if not text or text.startswith("skipped"):
        return []

    section_re = re.compile(
        r"=====\s*best_by_compleasm pass\s*\d+\s*\(([^)]+)\)\s*=====", re.IGNORECASE
    )
    sections = []
    matches = list(section_re.finditer(text))
    if matches:
        for i, m in enumerate(matches):
            label = m.group(1).strip()
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
            sections.append((label, text[start:end]))
    else:
        sections = [("", text)]

    def _num(pat, body):
        m = re.search(pat, body)
        if not m:
            return None
        try:
            return float(m.group(1))
        except ValueError:
            return None

    results = []
    for label, body in sections:
        braker_missing = _num(r"BRAKER is missing\s+([\d.]+)\s+BUSCOs", body)
        augustus_missing = _num(r"Augustus is missing\s+([\d.]+)\s+BUSCOs", body)
        genemark_missing = _num(r"GeneMark is missing\s+([\d.]+)\s+BUSCOs", body)

        decision = None
        if re.search(r"The BRAKER gene set .*? is the best one", body):
            decision = "kept the original merged gene set (already had the fewest missing BUSCOs)"
        elif re.search(r"All BUSCOs present in augustus.*?will be added to the braker", body):
            decision = "added BUSCO-supporting genes from AUGUSTUS and GeneMark on top of the merged set"
        elif re.search(r"Will enforce augustus.hints.gtf", body):
            decision = "rebuilt the gene set enforcing AUGUSTUS plus rescued BUSCOs"
        elif re.search(r"Will enforce genemark.gtf", body):
            decision = "rebuilt the gene set enforcing GeneMark plus rescued BUSCOs"
        elif re.search(r"there are no BUSCOs to be added", body):
            decision = "kept the original merged gene set (no rescuable BUSCOs)"
        elif re.search(r"WARNING: The new BRAKER gene set is not better", body):
            decision = "rebuilt the gene set, but it was not better than the original"

        new_missing = _num(r"It is missing\s+([\d.]+)% BUSCOs", body)

        if (braker_missing is None and augustus_missing is None and
                genemark_missing is None and decision is None):
            continue

        results.append({
            "label": label,
            "braker_missing": braker_missing,
            "augustus_missing": augustus_missing,
            "genemark_missing": genemark_missing,
            "decision": decision,
            "new_missing": new_missing,
        })
    return results


def format_bbc_decisions(passes, mode):
    """Render parsed best_by_compleasm passes as a methods-text sentence."""
    if not passes:
        return ""

    def _fmt(v):
        return f"{v:g}%" if v is not None else "n/a"

    def _one(p):
        bits = []
        bits.append(
            f"BRAKER {_fmt(p['braker_missing'])}, "
            f"AUGUSTUS {_fmt(p['augustus_missing'])}, "
            f"GeneMark {_fmt(p['genemark_missing'])} missing BUSCOs"
        )
        if p["decision"]:
            bits.append(p["decision"])
        if p["new_missing"] is not None:
            bits.append(f"final missing BUSCOs: {_fmt(p['new_missing'])}")
        return "; ".join(bits)

    if len(passes) == 1 and not passes[0]["label"]:
        return (
            "best_by_compleasm (Brůna, Gabriel & Hoff, 2025) compared the merged BRAKER, "
            "AUGUSTUS, and GeneMark gene sets using compleasm (Huang & Li, 2023) in protein "
            "mode and " + _one(passes[0]) + "."
        )

    pieces = []
    for p in passes:
        label = p["label"] or "pass"
        pieces.append(f"**{label}** — " + _one(p))
    return (
        "best_by_compleasm (Brůna, Gabriel & Hoff, 2025) was applied in three passes "
        "using compleasm (Huang & Li, 2023) in protein mode: " + " | ".join(pieces) + "."
    )


def generate_methods_text(workdir, mode):
    """Generate a single cohesive methods paragraph from pipeline outputs.

    Instead of collecting per-rule text chunks, this function inspects
    which output files exist and composes one methods section that reads
    like a paper's methods section with inline citations.
    """
    d = workdir
    parts = []

    def _count_unique_genes(gtf_path):
        """Count unique gene_id values — robust across GTF variants."""
        if not os.path.exists(gtf_path):
            return 0
        gids = set()
        with open(gtf_path) as f:
            for l in f:
                if l.startswith("#") or not l.strip():
                    continue
                m = re.search(r'gene_id "([^"]+)"', l)
                if m:
                    gids.add(m.group(1))
        return len(gids)

    # --- Pipeline + mode-specific citation ---
    pipeline_text = (
        "Gene prediction was performed using the BRAKER4 pipeline "
        "(Krall & Hoff, 2026) implemented in Snakemake (Mölder et al., 2021)"
    )
    if "BRAKER1" in mode:
        pipeline_text += ", running in ET/BRAKER1 mode (Hoff et al., 2016)."
    elif "BRAKER2" in mode:
        pipeline_text += ", running in EP/BRAKER2 mode (Bruna et al., 2021)."
    elif "BRAKER3" in mode:
        pipeline_text += ", running in ETP/BRAKER3 mode (Gabriel et al., 2024)."
        if "IsoSeq" in mode or "Dual" in mode:
            pipeline_text += " For IsoSeq integration, see also Bruna et al. (2025)."
    else:
        pipeline_text += "."
    parts.append(pipeline_text)

    # --- Repeat masking ---
    if os.path.exists(os.path.join(d, "preprocessing", "genome.fa.masked")):
        # Detect which masking tool was used from software_versions.tsv
        _versions_path = os.path.join(d, "software_versions.tsv")
        _used_red = False
        if os.path.exists(_versions_path):
            with open(_versions_path) as _vf:
                for _line in _vf:
                    if _line.strip().startswith("Red"):
                        _used_red = True
                        break
        if _used_red:
            parts.append(
                "The genome assembly was soft-masked for repeats using Red "
                "(Girgis, 2015), a neural-network-based repeat detector."
            )
        else:
            parts.append(
                "The genome assembly was soft-masked for repeats using RepeatModeler2 "
                "(Flynn et al., 2020) and RepeatMasker (Smit et al., 2013). "
                "Tandem repeats were additionally masked using Tandem Repeats Finder (Benson, 1999)."
            )

    # --- Alignment (must come before hint extraction) ---
    # Use software_versions.tsv as primary detection (survives cleanup)
    versions_file = os.path.join(d, "software_versions.tsv")
    recorded_tools = set()
    if os.path.exists(versions_file):
        for line in open(versions_file):
            tool = line.split('\t')[0].strip()
            recorded_tools.add(tool)

    # HISAT2 alignment
    has_hisat2 = "HISAT2" in recorded_tools or os.path.isdir(os.path.join(d, "hisat2_aligned"))
    if has_hisat2:
        parts.append(
            "Raw RNA-Seq reads were aligned to the genome using HISAT2 (Kim et al., 2019) "
            "with the --dta flag for downstream transcript assembly compatibility. "
            "Alignments were sorted and indexed with SAMtools (Danecek et al., 2021)."
        )

    # BAM preprocessing (when user provides pre-aligned BAMs without HISAT2)
    has_bam_sorted = os.path.isdir(os.path.join(d, "bam_sorted"))
    if has_bam_sorted and not has_hisat2:
        parts.append(
            "Pre-aligned RNA-Seq BAM files were validated, coordinate-sorted, "
            "and indexed using SAMtools (Danecek et al., 2021)."
        )

    # IsoSeq alignment
    has_minimap2 = "minimap2" in recorded_tools or os.path.isdir(os.path.join(d, "minimap2_aligned")) or os.path.isdir(os.path.join(d, "isoseq_aligned"))
    has_isoseq_sorted = os.path.isdir(os.path.join(d, "isoseq_sorted")) or os.path.isdir(os.path.join(d, "isoseq_merged"))
    if has_minimap2:
        parts.append(
            "PacBio IsoSeq long reads were aligned to the genome using minimap2 "
            "(Li, 2018) in splice:hq mode. "
            "Alignments were sorted and indexed with SAMtools (Danecek et al., 2021)."
        )
    elif has_isoseq_sorted:
        parts.append(
            "Pre-aligned PacBio IsoSeq BAM files were sorted and indexed "
            "using SAMtools (Danecek et al., 2021)."
        )

    # --- RNA-Seq hint extraction (after alignment) ---
    hints_file = os.path.join(d, "bam2hints.gff")
    if os.path.exists(hints_file):
        n_hints = sum(1 for _ in open(hints_file))
        parts.append(
            f"Spliced alignments were converted to {n_hints:,} intron hints "
            "using bam2hints from AUGUSTUS (Stanke et al., 2008). "
            "Strand orientation was determined and spurious introns filtered "
            "using filterIntronsFindStrand.pl based on canonical splice site dinucleotides."
        )

    # --- VARUS ---
    varus_stats = os.path.join(d, "varus", "varus_stats.txt")
    varus_runlist = os.path.join(d, "varus", "varus_runlist.tsv")
    if os.path.isdir(os.path.join(d, "varus")):
        varus_text = (
            "RNA-Seq data was automatically selected and downloaded from the NCBI "
            "Sequence Read Archive using VARUS (Stanke et al., 2019)"
        )
        if os.path.exists(varus_stats):
            stats_text = read_file(varus_stats)
            m = re.search(r"SRA experiments downloaded:\s*(\d+)", stats_text)
            if m:
                varus_text += f", which selected **{int(m.group(1)):,} SRA experiments**"
            m = re.search(r"Mapped reads:\s*(\d+)", stats_text)
            if m:
                varus_text += f" yielding **{int(m.group(1)):,} mapped reads**"
        varus_text += "."
        parts.append(varus_text)

    # --- SRA download ---
    if os.path.isdir(os.path.join(d, "sra_fastq")):
        parts.append(
            "RNA-Seq data was downloaded from the NCBI Sequence Read Archive "
            "using the SRA Toolkit (SRA Toolkit Development Team, 2020)."
        )

    # --- ProtHint ---
    if os.path.exists(os.path.join(d, "prothint_hints.gff")):
        n_phints = sum(1 for _ in open(os.path.join(d, "prothint_hints.gff")))
        parts.append(
            f"Protein evidence was generated using ProtHint (Bruna et al., 2020), "
            f"which produced {n_phints:,} protein-based hints by aligning the protein "
            "database to the genome using DIAMOND (Buchfink et al., 2021) and "
            "Spaln (Iwata & Gotoh, 2012)."
        )
    if os.path.exists(os.path.join(d, "prothint_hints_iter2.gff")):
        parts.append(
            "ProtHint was re-run in a second iteration using the initial AUGUSTUS "
            "predictions as additional seeds, producing refined protein hints."
        )

    # --- GeneMark ---
    if "ES" in mode and "ab initio" in mode:
        n = _count_unique_genes(os.path.join(d, "GeneMark-ES", "genemark.gtf"))
        if n > 0:
            parts.append(
                f"Ab initio gene prediction was performed using GeneMark-ES "
                f"(Lomsadze et al., 2005), which predicted {n:,} genes."
            )
    elif "ET" in mode and "BRAKER1" in mode:
        n = _count_unique_genes(os.path.join(d, "genemark", "genemark.gtf"))
        if n > 0:
            parts.append(
                f"GeneMark-ET (Lomsadze et al., 2014) was trained on RNA-Seq intron "
                f"hints and predicted {n:,} genes."
            )
    elif "EP" in mode and "BRAKER2" in mode:
        n = _count_unique_genes(os.path.join(d, "GeneMark-EP", "genemark.gtf"))
        if n > 0:
            parts.append(
                f"GeneMark-EP+ (Bruna et al., 2020) was trained with protein-derived "
                f"intron hints and predicted {n:,} genes."
            )
    elif "BRAKER3" in mode:
        gm_sr_gtf = os.path.join(d, "GeneMark-ETP", "genemark.gtf")
        gm_iso_gtf = os.path.join(d, "GeneMark-ETP-isoseq", "genemark.gtf")
        is_dual_etp = os.path.exists(gm_sr_gtf) and os.path.exists(gm_iso_gtf)

        if is_dual_etp:
            n_sr = _count_unique_genes(gm_sr_gtf)
            n_iso = _count_unique_genes(gm_iso_gtf)
            parts.append(
                f"In dual mode, GeneMark-ETP (Bruna et al., 2024) was run twice in parallel: "
                f"once on the short-read RNA-Seq + protein evidence "
                f"(predicting **{n_sr:,} genes**) and once on the IsoSeq long-read + protein "
                f"evidence (predicting **{n_iso:,} genes**). Both runs internally employed "
                f"GeneMarkS-T (Tang et al., 2015) for transcript-based gene finding and "
                f"DIAMOND (Buchfink et al., 2021) for protein alignment."
            )
            n_merged = _count_unique_genes(os.path.join(d, "dual_etp_merged", "training.gtf"))
            if n_merged > 0:
                parts.append(
                    f"The high-confidence training genes from both GeneMark-ETP runs were "
                    f"merged and deduplicated with TSEBRA (Gabriel et al., 2021), yielding "
                    f"**{n_merged:,} training genes** for AUGUSTUS."
                )
        else:
            n = _count_unique_genes(gm_sr_gtf)
            if n > 0:
                # The evidence-source wording depends on whether this is
                # IsoSeq-only + proteins or short-read RNA-Seq + proteins.
                # detect_mode() already distinguishes the two.
                if "IsoSeq" in mode:
                    evidence_desc = "PacBio IsoSeq long reads and protein evidence"
                else:
                    evidence_desc = "short-read RNA-Seq and protein evidence"
                parts.append(
                    f"GeneMark-ETP (Bruna et al., 2024) was trained using {evidence_desc}, "
                    f"internally employing GeneMarkS-T (Tang et al., 2015) "
                    f"for transcript-based gene finding and DIAMOND (Buchfink et al., 2021) "
                    f"for protein alignment, and predicted {n:,} genes."
                )

    # --- Training gene selection pipeline ---
    # GeneMark filtering step. The source file and phrasing depend on mode:
    #   - ES: no filtering happens (genemark.f.good.gtf is a symlink to the
    #         raw predictions); skip the paragraph entirely.
    #   - ET/EP: filter_genemark.smk runs filterGenemark.pl against our hints
    #         and writes genemark/genemark.f.good.gtf (or GeneMark-EP/...).
    #         Describe it as our filtering step.
    #   - ETP/IsoSeq/dual: GeneMark-ETP itself produces a pre-filtered
    #         high-confidence training gene set as training.gtf. Our own
    #         filter_genemark rule is NOT in the DAG for these modes, so
    #         genemark/genemark.f.good.gtf does not exist. Read the ETP
    #         training.gtf directly and describe it as ETP-internal filtering.
    if "ES" in mode and "ab initio" in mode:
        pass  # ES has no evidence to filter by
    elif "Dual" in mode:
        # Dual mode already announced the merged training-gene count in the
        # TSEBRA-merge paragraph above (dual_etp_merged/training.gtf). Adding
        # an "ETP-internal high-confidence selection" line here would repeat
        # the same number and read as duplicated narrative.
        pass
    elif "ETP" in mode or "IsoSeq" in mode:
        # Single-run ETP (short-read or IsoSeq): describe the ETP-internal
        # high-confidence training gene selection from GeneMark-ETP/training.gtf.
        p = os.path.join(d, "GeneMark-ETP", "training.gtf")
        if os.path.exists(p):
            n_good = _count_unique_genes(p)
            parts.append(
                f"GeneMark-ETP's internal high-confidence training gene selection "
                f"yielded {n_good:,} training gene candidates."
            )
    else:
        # ET / EP — our own filter_genemark rule wrote these paths
        for gff_path in [
            os.path.join(d, "genemark", "genemark.f.good.gtf"),
            os.path.join(d, "GeneMark-EP", "genemark.f.good.gtf"),
        ]:
            if os.path.exists(gff_path):
                n_good = _count_unique_genes(gff_path)
                parts.append(
                    f"GeneMark predictions were filtered for genes with evidence support, "
                    f"yielding {n_good:,} high-quality training gene candidates."
                )
                break

    # Etraining filtering — runs BEFORE DIAMOND redundancy removal (post
    # BRAKER_PARITY_AUDIT.md Cluster A fix #7). The report must describe the
    # pipeline steps in their actual execution order: etraining → DIAMOND.
    etrain_file = os.path.join(d, "etraining_gene_count.txt")
    if os.path.exists(etrain_file):
        text = read_file(etrain_file).strip()
        m = re.match(r"(\d+)\s*->\s*(\d+)", text)
        if m and m.group(1) != m.group(2):
            parts.append(
                f"After etraining filtering, {int(m.group(2)):,} of "
                f"{int(m.group(1)):,} genes were retained."
            )

    # Redundancy removal (runs AFTER etraining filtering)
    locus_file = os.path.join(d, "locus_count.txt")
    if os.path.exists(locus_file):
        text = read_file(locus_file).strip()
        m = re.match(r"(\d+)\s*->\s*(\d+)", text)
        if m:
            parts.append(
                f"Redundant sequences were then removed using DIAMOND (Buchfink et al., 2021), "
                f"reducing the training set from {int(m.group(1)):,} to "
                f"{int(m.group(2)):,} loci."
            )

    # Poisson downsampling (runs AFTER DIAMOND). Without this paragraph the
    # report jumps straight from "reduced to N loci" to the train/test split,
    # leaving the reader to wonder why N doesn't equal train+test when the
    # downsampler dropped genes in between.
    #
    # The downsample rule writes one of two formats to the counter file:
    #   "N -> M (removed: ...)"      — normal case
    #   "Skipped: X% single-exon genes" — 80%-single-exon skip path
    # We handle both so the reader always sees what happened.
    ds_file = os.path.join(d, "downsample_gene_count.txt")
    if os.path.exists(ds_file):
        text = read_file(ds_file).strip()
        m = re.match(r"(\d+)\s*->\s*(\d+)", text)
        m_skip = re.match(r"Skipped:\s*(\d+)%\s*single-exon", text)
        if m and m.group(1) != m.group(2):
            parts.append(
                f"Training genes were then downsampled using a Poisson distribution "
                f"(downsample_traingenes.pl, lambda=2) to reduce single-exon gene "
                f"over-representation, yielding {int(m.group(2)):,} of "
                f"{int(m.group(1)):,} genes."
            )
        elif m_skip:
            parts.append(
                f"Poisson downsampling was skipped because {m_skip.group(1)}% of the "
                f"training genes are single-exon (BRAKER4 heuristic: skip downsampling "
                f"when the majority of genes are single-exon)."
            )

    # Train/test split
    split_log = os.path.join(d, "training_set_split.log")
    if os.path.exists(split_log):
        text = read_file(split_log)
        m_train = re.search(r"Training set.*?(\d+)\s*genes", text)
        m_test = re.search(r"Test set.*?(\d+)\s*genes", text)
        if m_train and m_test:
            parts.append(
                f"The training set was split into {int(m_train.group(1)):,} genes "
                f"for training and {int(m_test.group(1)):,} genes held out for "
                f"accuracy assessment."
            )

    # --- AUGUSTUS training ---
    acc_before = os.path.join(d, "accuracy_after_training.txt")
    acc_after = os.path.join(d, "accuracy_after_optimize.txt")
    if os.path.exists(acc_before):
        optimized = False
        if os.path.exists(acc_after):
            acc_text = read_file(acc_after)
            if "SKIPPED" not in acc_text and "skipped" not in acc_text:
                optimized = True
        if optimized:
            parts.append(
                "AUGUSTUS (Stanke et al., 2008) species-specific parameters were trained "
                "on the selected gene set using etraining and further refined "
                "using optimize_augustus.pl."
            )
        else:
            parts.append(
                "AUGUSTUS (Stanke et al., 2008) species-specific parameters were trained "
                "on the selected gene set using etraining."
            )

    # --- AUGUSTUS prediction ---
    aug_gtf = os.path.join(d, "augustus.hints.gtf")
    if os.path.exists(aug_gtf):
        n = sum(1 for l in open(aug_gtf) if "\tgene\t" in l)
        if "ES" in mode and "ab initio" in mode:
            parts.append(
                f"AUGUSTUS predicted {n:,} genes using the trained parameters "
                "with BUSCO-based hints from compleasm (Huang & Li, 2023) as the only extrinsic evidence."
            )
        else:
            parts.append(
                f"AUGUSTUS predicted {n:,} genes using the trained parameters and "
                "all available extrinsic evidence as hints."
            )

    # --- TSEBRA ---
    tsebra_gtf = os.path.join(d, "braker.tsebra.gtf")
    if os.path.exists(tsebra_gtf):
        parts.append(
            "Gene predictions from AUGUSTUS and GeneMark were merged using TSEBRA "
            "(Gabriel et al., 2021)."
        )

    # --- best_by_compleasm decisions (parsed from the aggregated log) ---
    bbc_log_path = os.path.join(d, "best_by_compleasm.log")
    bbc_passes = parse_bbc_log(bbc_log_path)
    if bbc_passes:
        parts.append(format_bbc_decisions(bbc_passes, mode))

    # --- Final gene set ---
    braker_gtf = os.path.join(d, "braker.gtf")
    if os.path.exists(braker_gtf):
        n_genes = sum(1 for l in open(braker_gtf) if "\tgene\t" in l)
        n_tx = sum(1 for l in open(braker_gtf) if "\ttranscript\t" in l or "\tmRNA\t" in l)
        parts.append(
            f"The final gene set contains **{n_genes:,} genes** and "
            f"**{n_tx:,} transcripts** after filtering for internal stop codons "
            "and CDS normalization."
        )

    # --- UTR ---
    utr_gtf = os.path.join(d, "braker_utr.gtf")
    def _is_utr_line(line):
        # AUGUSTUS uses "5'-UTR"/"3'-UTR", stringtie2utr.py uses "five_prime_UTR"/"three_prime_UTR"
        low = line.lower()
        return ("\t5'-utr\t" in low or "\tfive_prime_utr\t" in low or
                "\t3'-utr\t" in low or "\tthree_prime_utr\t" in low)

    if os.path.exists(utr_gtf):
        n_utr_genes = sum(1 for l in open(utr_gtf) if _is_utr_line(l))
        # Count unique gene IDs with UTR features
        utr_gene_ids = set()
        for l in open(utr_gtf):
            if _is_utr_line(l):
                import re as _re
                m = _re.search(r'gene_id "([^"]+)"', l)
                if m:
                    utr_gene_ids.add(m.group(1))
        n_genes_with_utr = len(utr_gene_ids)
        if n_genes_with_utr > 0:
            parts.append(
                f"UTR features were added by decorating gene models with StringTie2 "
                f"(Kovaka et al., 2019) transcript assemblies. "
                f"**{n_genes_with_utr:,}** genes received UTR annotations "
                f"({n_utr_genes:,} UTR features total)."
            )

    # --- rRNA ---
    rrna_gff = os.path.join(d, "barrnap", "rRNA.gff3")
    if os.path.exists(rrna_gff):
        n_rrna = sum(1 for l in open(rrna_gff) if not l.startswith("#") and l.strip())
        if n_rrna > 0:
            parts.append(
                f"**{n_rrna:,}** ribosomal RNA genes were predicted using barrnap (Seemann, 2018) "
                "and merged into the GFF3 annotation."
            )
        else:
            parts.append(
                "No ribosomal RNA genes were detected by barrnap (Seemann, 2018) in this genome."
            )

    # --- ncRNA (tRNAscan-SE + Infernal) ---
    ncrna_parts = []
    trna_gff = os.path.join(d, "ncrna", "tRNAs.gff3")
    if os.path.exists(trna_gff):
        n_trna = sum(1 for l in open(trna_gff) if not l.startswith("#") and l.strip())
        if n_trna > 0:
            ncrna_parts.append(
                f"**{n_trna:,}** transfer RNA genes were predicted using tRNAscan-SE "
                "(Chan & Lowe, 2019)."
            )

    infernal_gff = os.path.join(d, "ncrna", "ncRNAs_infernal.gff3")
    if os.path.exists(infernal_gff):
        n_inf = sum(1 for l in open(infernal_gff) if not l.startswith("#") and l.strip())
        if n_inf > 0:
            ncrna_parts.append(
                f"**{n_inf:,}** non-coding RNA features (snoRNAs, snRNAs, miRNAs, etc.) "
                "were identified by scanning against the Rfam database "
                "(Kalvari et al., 2021) using Infernal (Nawrocki & Eddy, 2013)."
            )

    lncrna_gff = os.path.join(d, "ncrna", "lncRNAs.gff3")
    if os.path.exists(lncrna_gff):
        n_lnc = sum(1 for l in open(lncrna_gff) if not l.startswith("#") and l.strip())
        if n_lnc > 0:
            ncrna_parts.append(
                f"**{n_lnc:,}** long non-coding RNAs were identified from StringTie "
                "transcriptome assemblies using FEELnc (Wucher et al., 2017)."
            )

    if ncrna_parts:
        ncrna_parts.append(
            "All ncRNA predictions were merged into the final GFF3 annotation."
        )
        parts.extend(ncrna_parts)

    # --- GFF3 ---
    if os.path.exists(os.path.join(d, "braker.gff3")):
        parts.append(
            "The final annotation was converted to GFF3 format using AGAT "
            "(Dainat, 2024)."
        )

    # --- FANTASIA ---
    fantasia_summary = os.path.join(d, "fantasia", "fantasia_summary.txt")
    if os.path.exists(fantasia_summary):
        n_proteins_hc = None
        try:
            saw_proteins_line = False
            with open(fantasia_summary) as f:
                for line in f:
                    if line.lstrip().startswith("Proteins with at least one GO term"):
                        saw_proteins_line = True
                        continue
                    if saw_proteins_line and "above cutoff" in line:
                        n_proteins_hc = line.split(":")[-1].strip().replace(",", "")
                        break
        except Exception:
            n_proteins_hc = None
        if n_proteins_hc and n_proteins_hc.isdigit():
            parts.append(
                f"Predicted proteins were functionally annotated with FANTASIA-Lite "
                f"(Mart\u00ednez-Redondo et al., 2025; Cases et al., 2025), which assigns "
                f"GO terms via ProtT5 (Rostlab/prot_t5_xl_uniref50) protein language model "
                f"embeddings against a bundled lookup of pre-computed reference embeddings; "
                f"**{int(n_proteins_hc):,}** proteins received at least one high-confidence "
                f"GO assignment. The GO terms were merged into the GFF3 annotation as the "
                f"reserved <code>Ontology_term</code> attribute on each mRNA feature, with "
                f"the per-gene union rolled up to the parent gene feature."
            )
        else:
            parts.append(
                "Predicted proteins were functionally annotated with FANTASIA-Lite "
                "(Mart\u00ednez-Redondo et al., 2025; Cases et al., 2025), which assigns "
                "GO terms via ProtT5 protein language model embeddings against a bundled "
                "lookup of pre-computed reference embeddings. The GO terms were merged "
                "into the GFF3 annotation as the reserved <code>Ontology_term</code> "
                "attribute on each mRNA feature, with the per-gene union rolled up to the "
                "parent gene feature."
            )

    # --- QC ---
    qc_parts = []
    if os.path.isdir(os.path.join(d, "busco")):
        qc_parts.append("BUSCO (Manni et al., 2021)")
    if os.path.isdir(os.path.join(d, "compleasm_proteins")) or os.path.isdir(os.path.join(d, "compleasm_genome_out")):
        qc_parts.append("compleasm (Huang & Li, 2023) with miniprot (Li, 2023)")
    if os.path.isdir(os.path.join(d, "omark")):
        qc_parts.append("OMArk (Nevers et al., 2025) with OMAmer (Rossier et al., 2021)")
    if os.path.isdir(os.path.join(d, "gffcompare")):
        qc_parts.append("gffcompare (Pertea & Pertea, 2020)")
    if qc_parts:
        parts.append(
            "Quality was assessed using " + ", ".join(qc_parts) + "."
        )

    return " ".join(parts)


def deduplicate_citations(text):
    """Deduplicate plain-text citations (one per line)."""
    seen = set()
    result = []
    for line in text.strip().splitlines():
        line = line.strip()
        if line and line not in seen:
            seen.add(line)
            result.append(line)
    return result


def deduplicate_bibtex(text):
    """Deduplicate BibTeX entries by key."""
    entries = OrderedDict()
    current_key = None
    current_entry = []

    for line in text.splitlines():
        m = re.match(r"@\w+\{(\S+),", line)
        if m:
            if current_key and current_entry:
                entries[current_key] = "\n".join(current_entry)
            current_key = m.group(1)
            current_entry = [line]
        else:
            current_entry.append(line)

    if current_key and current_entry:
        entries[current_key] = "\n".join(current_entry)

    return "\n\n".join(entries.values())


def format_compleasm_as_busco(text):
    """Convert compleasm multi-line summary to BUSCO-style one-liner.

    Input format:
        ## lineage: eukaryota_odb12
        S:1.55%, 2
        D:0.00%, 0
        F:0.78%, 1
        I:0.00%, 0
        M:97.67%, 126
        N:129

    Output format:
        Lineage: eukaryota_odb12
        C:1.6%[S:1.6%,D:0.0%],F:0.8%,I:0.0%,M:97.7%,n:129
    """
    vals = {}
    lineage = ""
    for line in text.strip().splitlines():
        line = line.strip()
        m = re.match(r"##\s*lineage:\s*(\S+)", line)
        if m:
            lineage = m.group(1)
            continue
        m = re.match(r"([SDFIMN]):?([\d.]+)%?,?\s*(\d+)?", line)
        if m:
            key = m.group(1)
            if key == "N":
                vals["N"] = m.group(2)
            else:
                vals[key] = float(m.group(2))
    if "S" not in vals or "D" not in vals:
        return text  # can't parse, return as-is
    c = vals.get("S", 0) + vals.get("D", 0)
    parts = [f"Lineage: {lineage}"] if lineage else []
    parts.append(
        f"\tC:{c:.1f}%[S:{vals.get('S', 0):.1f}%,D:{vals.get('D', 0):.1f}%],"
        f"F:{vals.get('F', 0):.1f}%,I:{vals.get('I', 0):.1f}%,"
        f"M:{vals.get('M', 0):.1f}%,n:{vals.get('N', '?')}"
    )
    return "\n".join(parts)


def read_qc_summary(path, reformat_compleasm=False):
    """Read a QC summary file and return as HTML paragraph."""
    text = read_file(path).strip()
    if not text:
        return ""
    if reformat_compleasm:
        text = format_compleasm_as_busco(text)
    return f"<pre>{text}</pre>"


def collect_benchmarks(workdir):
    """Collect Snakemake benchmark files from benchmarks/{sample}/.

    Returns list of dicts with keys: rule, wallclock_s, cpu_s, max_rss_mb.
    Snakemake benchmark format (TSV, 1 header + 1 data row):
      s  h:m:s  max_rss  max_vms  max_uss  max_pss  io_in  io_out  mean_load  cpu_time
    """
    import glob as globmod
    benchmarks = []

    # Find benchmark dir — could be benchmarks/{sample}/ or ../benchmarks/{sample}/
    bench_dirs = []
    parent = os.path.dirname(workdir.rstrip("/"))
    sample = os.path.basename(workdir.rstrip("/"))
    for candidate in [
        os.path.join(parent, "benchmarks", sample),
        os.path.join(workdir, "..", "benchmarks", sample),
        os.path.join(os.path.dirname(parent), "benchmarks", sample),
    ]:
        if os.path.isdir(candidate):
            bench_dirs.append(candidate)
            break

    if not bench_dirs:
        return benchmarks

    bench_dir = bench_dirs[0]
    for txt in sorted(globmod.glob(os.path.join(bench_dir, "*", "*.txt"))):
        dir_name = os.path.basename(os.path.dirname(txt))
        file_stem = os.path.splitext(os.path.basename(txt))[0]
        # If the filename differs from the directory name, use the filename
        # (this keeps BUSCO genome vs proteome distinguishable)
        rule_name = file_stem if file_stem != dir_name else dir_name
        try:
            with open(txt) as f:
                header = f.readline().strip().split("\t")
                values = f.readline().strip().split("\t")
            if len(values) < 2:
                continue
            row = dict(zip(header, values))
            wallclock = float(row.get("s", 0))
            cpu_time = float(row.get("cpu_time", 0))
            max_rss = float(row.get("max_rss", 0))  # in MB
            benchmarks.append({
                "rule": rule_name,
                "wallclock_s": wallclock,
                "cpu_s": cpu_time,
                "max_rss_mb": max_rss,
            })
        except (ValueError, IndexError, KeyError):
            continue

    # Sort by wallclock descending
    benchmarks.sort(key=lambda x: x["wallclock_s"], reverse=True)
    return benchmarks


def format_time(seconds):
    """Format seconds into human-readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}min"
    else:
        return f"{seconds/3600:.1f}h"


def generate_benchmark_plot(benchmarks, outdir):
    """Generate publication-quality resource consumption plots.

    Creates a two-panel figure:
      Left:  Wall clock vs CPU time (minutes) — horizontal bars
      Right: Peak RAM usage (GB) — horizontal bars

    Returns the path to the PNG file, or None if matplotlib is unavailable.
    """
    if not benchmarks:
        return None

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        return None

    # Filter to rules that took > 1 second
    data = [b for b in benchmarks if b["wallclock_s"] > 1]
    if not data:
        return None

    # Top 20 by wallclock
    data = data[:20]
    data.reverse()  # longest at top

    rules = [d["rule"].replace("_", " ") for d in data]
    wallclock = [d["wallclock_s"] / 60 for d in data]  # minutes
    cpu = [d["cpu_s"] / 60 for d in data]
    ram = [d["max_rss_mb"] / 1024 for d in data]  # GB

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(4, len(rules) * 0.38)),
                                    gridspec_kw={"width_ratios": [3, 2]})

    y = range(len(rules))
    bar_h = 0.38

    # Left panel: wall clock + CPU time
    ax1.barh(y, wallclock, color="#4C72B0", alpha=0.85, label="Wall clock",
             height=bar_h, align="edge")
    ax1.barh([i - bar_h for i in y], cpu, color="#DD8452", alpha=0.85,
             label="CPU time", height=bar_h, align="edge")
    ax1.set_yticks(y)
    ax1.set_yticklabels(rules, fontsize=8)
    ax1.set_xlabel("Time (minutes)", fontsize=10)
    ax1.set_title("Runtime", fontsize=12, fontweight="bold")
    ax1.legend(fontsize=8, loc="lower right", framealpha=0.9)
    ax1.grid(axis="x", alpha=0.3, linewidth=0.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Right panel: peak RAM
    ax2.barh(y, ram, color="#55A868", alpha=0.85, height=bar_h * 2)
    ax2.set_yticks(y)
    ax2.set_yticklabels([], fontsize=8)  # labels already on left panel
    ax2.set_xlabel("Peak RAM (GB)", fontsize=10)
    ax2.set_title("Memory", fontsize=12, fontweight="bold")
    ax2.grid(axis="x", alpha=0.3, linewidth=0.5)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

    plt.tight_layout()
    plot_path = os.path.join(outdir, "runtime_plot.png")
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    return plot_path


def detect_mode(workdir):
    """Detect BRAKER mode from directory structure.

    Order matters: check specific modes first (ETP > EP > ES),
    then fall back to ET. The genemark/ directory exists in all modes
    (from filter_genemark), so it must be checked last.

    Returns a label that includes the BRAKER version name:
    - ET  = BRAKER1 (Hoff et al., 2016)
    - EP  = BRAKER2 (Bruna et al., 2021)
    - ETP = BRAKER3 (Gabriel et al., 2024)
    - IsoSeq/Dual: cite the book chapter (Hoff et al., 2024)
    """
    if os.path.exists(os.path.join(workdir, "GeneMark-ETP")):
        if os.path.exists(os.path.join(workdir, "GeneMark-ETP-isoseq")):
            return "Dual/BRAKER3 (short-read + IsoSeq + proteins)"
        # Check for IsoSeq-only (no short-read bam2hints)
        if not os.path.exists(os.path.join(workdir, "bam2hints.gff")) and \
           not os.path.exists(os.path.join(workdir, "genemark")):
            return "IsoSeq/BRAKER3 (long-read + proteins)"
        return "ETP/BRAKER3 (RNA-Seq + proteins)"
    elif os.path.exists(os.path.join(workdir, "GeneMark-EP")):
        return "EP/BRAKER2 (proteins only)"
    elif os.path.exists(os.path.join(workdir, "GeneMark-ES")):
        if not os.path.exists(os.path.join(workdir, "bam2hints.gff")) and \
           not os.path.exists(os.path.join(workdir, "prothint_hints.gff")):
            return "ES (ab initio)"
        return "EP/BRAKER2 (proteins only)"
    elif os.path.exists(os.path.join(workdir, "genemark")):
        return "ET/BRAKER1 (RNA-Seq only)"
    return "Unknown"


CSS = """
body {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    max-width: 900px;
    margin: 40px auto;
    padding: 0 20px;
    color: #333;
    line-height: 1.6;
    background: #fafafa;
}
h1 {
    color: #2c5f2d;
    border-bottom: 3px solid #2c5f2d;
    padding-bottom: 10px;
}
h2 {
    color: #3a7d44;
    border-bottom: 1px solid #ccc;
    padding-bottom: 5px;
    margin-top: 30px;
}
h3 { color: #555; }
.meta {
    background: #e8f5e9;
    border-left: 4px solid #2c5f2d;
    padding: 12px 16px;
    margin: 20px 0;
    font-size: 0.95em;
}
.meta dt { font-weight: bold; display: inline; }
.meta dd { display: inline; margin-left: 5px; }
.meta dd::after { content: ''; display: block; }
pre {
    background: #f5f5f5;
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 12px;
    overflow-x: auto;
    font-size: 0.85em;
}
.citation-list {
    font-size: 0.9em;
    line-height: 1.8;
}
.citation-list li { margin-bottom: 6px; }
table {
    border-collapse: collapse;
    width: 100%;
    margin: 10px 0;
}
th, td {
    border: 1px solid #ddd;
    padding: 8px 12px;
    text-align: left;
}
th { background: #f0f0f0; }
.figure {
    text-align: center;
    margin: 20px 0;
}
.figure img { border: 1px solid #ddd; border-radius: 4px; }
.figure .caption {
    font-size: 0.85em;
    color: #666;
    margin-top: 6px;
    font-style: italic;
}
details { margin: 5px 0; }
summary {
    cursor: pointer;
    font-weight: bold;
    color: #3a7d44;
}
footer {
    margin-top: 40px;
    padding-top: 10px;
    border-top: 1px solid #ccc;
    font-size: 0.8em;
    color: #888;
}
"""


def generate_html(sample_name, mode, methods_text, citations, images, qc_data, benchmarks=None, outdir=".", software_versions=None):
    """Generate the full HTML report."""

    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Embed logo if available
    logo_html = ""
    for logo_candidate in [
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "img", "logo.jpeg"),
        os.path.join(outdir, "..", "..", "img", "logo.jpeg"),
    ]:
        if os.path.exists(logo_candidate):
            logo_html = f'<div style="text-align:center; margin-bottom:10px;">{embed_image(logo_candidate)}</div>'
            # Make logo smaller
            logo_html = logo_html.replace('max-width:100%', 'max-width:120px')
            break

    # Pipeline version (from software_versions.tsv)
    pipeline_version = ""
    if software_versions:
        for tool, ver in software_versions:
            if tool == "BRAKER4 pipeline":
                pipeline_version = ver
                break

    # Genome location
    genome_location = ""
    genome_path_file = os.path.join(outdir, "genome_path.txt")
    if os.path.exists(genome_path_file):
        genome_location = read_file(genome_path_file).strip()

    # Build result file listing for header
    # Files will be gzipped AFTER report generation, so check both .gz and uncompressed
    result_files_html = ""
    result_files = [
        ("braker.gtf", "Gene predictions (GTF)"),
        ("braker.gff3", "Gene predictions (GFF3)"),
        ("braker_with_ncRNA.gff3", "Genes + ncRNA (GFF3; rRNA, tRNA, snoRNA, snRNA, miRNA, lncRNA — only when run_ncrna=1)"),
        ("braker.go.gff3", "Gene predictions with FANTASIA GO terms (GFF3 Ontology_term attribute on mRNA + gene features — only when run_fantasia=1)"),
        ("braker_with_ncRNA.go.gff3", "Genes + ncRNA + FANTASIA GO terms (GFF3 — only when run_fantasia=1 and run_ncrna=1)"),
        ("fantasia_go_terms.tsv", "FANTASIA per-transcript GO term table (transcript_id, go_id, go_name, go_namespace, reliability_index — only when run_fantasia=1)"),
        ("braker.aa", "Protein sequences (FASTA)"),
        ("braker.codingseq", "Coding sequences (FASTA)"),
        ("braker_utr.gtf", "Gene predictions with UTR features"),
        ("genome.fa", "Repeat-masked genome assembly (pipeline-generated)"),
        ("hintsfile.gff", "All extrinsic evidence hints"),
        ("gene_support.tsv", "Per-gene evidence support statistics"),
        ("varus_runlist.tsv", "VARUS: SRA experiments selected and downloaded"),
        ("varus_stats.txt", "VARUS: Run summary statistics"),
        ("braker_report.html", "This report"),
        ("braker_citations.bib", "BibTeX citations for all tools used"),
    ]
    for filename, description in result_files:
        # Check for .gz version first, then uncompressed
        gz_path = os.path.join(outdir, filename + ".gz")
        plain_path = os.path.join(outdir, filename)
        if os.path.exists(gz_path):
            filepath = gz_path
            display_name = filename + ".gz"
        elif os.path.exists(plain_path):
            filepath = plain_path
            # These will be gzipped after report generation
            display_name = filename + ".gz" if filename not in ("gene_support.tsv", "braker_report.html", "braker_citations.bib") else filename
        else:
            continue

        # Skip braker_utr.gtf if no UTR features were actually added
        if filename == "braker_utr.gtf":
            has_utr = False
            try:
                if filepath.endswith(".gz"):
                    import gzip
                    opener = lambda p: gzip.open(p, "rt")
                else:
                    opener = lambda p: open(p)
                with opener(filepath) as fh:
                    for line in fh:
                        low = line.lower()
                        if ("\t5'-utr\t" in low or "\tfive_prime_utr\t" in low or
                            "\t3'-utr\t" in low or "\tthree_prime_utr\t" in low):
                            has_utr = True
                            break
            except Exception:
                has_utr = True  # if check fails, keep the file listed
            if not has_utr:
                continue
        size_bytes = os.path.getsize(filepath)
        if size_bytes > 1024 * 1024:
            size_str = f"{size_bytes / (1024*1024):.1f} MB"
        elif size_bytes > 1024:
            size_str = f"{size_bytes / 1024:.0f} KB"
        else:
            size_str = f"{size_bytes} B"
        result_files_html += f'<tr><td><a href="{display_name}"><code>{display_name}</code></a> ({size_str})</td><td>{description}</td></tr>\n'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>BRAKER4 Report: {sample_name}</title>
<style>{CSS}</style>
</head>
<body>
{logo_html}
<h1>BRAKER4 Gene Prediction Report</h1>

<div class="meta">
<dl>
<dt>Pipeline version:</dt><dd>BRAKER4 v{pipeline_version}</dd>
<dt>Sample:</dt><dd>{sample_name}</dd>
<dt>Mode:</dt><dd>{mode}</dd>
<dt>Genome:</dt><dd><code>{genome_location}</code></dd>
<dt>Date:</dt><dd>{now}</dd>
</dl>
</div>

<h2>Output Files</h2>
<table>
<tr><th>File</th><th>Description</th></tr>
{result_files_html}
</table>
"""

    # Methods section — one cohesive paragraph
    if methods_text:
        methods_html = methods_text
        methods_html = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", methods_html)
        methods_html = re.sub(r"\*(.+?)\*", r"<em>\1</em>", methods_html)
        html += f"\n<h2>Methods</h2>\n<p>{methods_html}</p>\n"

    # VARUS run list (if available) — collapsible table of SRA experiments
    varus_runlist_path = None
    for search_dir in [outdir, os.path.join(outdir, "..")]:
        p = os.path.join(search_dir, "varus_runlist.tsv")
        if os.path.exists(p):
            varus_runlist_path = p
            break
    if varus_runlist_path:
        html += "\n<h2>VARUS: Selected SRA Experiments</h2>\n"
        html += "<p>The following SRA experiments were automatically selected by VARUS:</p>\n"
        html += "<details open>\n<summary>SRA run list</summary>\n"
        html += "<table>\n<tr><th>SRA Accession</th><th>Reads Downloaded</th><th>Status</th></tr>\n"
        with open(varus_runlist_path) as f:
            header = f.readline()  # skip header
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    accession = cols[0]
                    reads = cols[1]
                    status = cols[2]
                    # Make accession a link to NCBI
                    link = f'<a href="https://www.ncbi.nlm.nih.gov/sra/{accession}" target="_blank">{accession}</a>'
                    html += f"<tr><td>{link}</td><td>{int(reads):,}</td><td>{status}</td></tr>\n"
        html += "</table>\n</details>\n"

        # VARUS stats summary
        varus_stats_path = os.path.join(os.path.dirname(varus_runlist_path), "varus_stats.txt")
        if os.path.exists(varus_stats_path):
            html += read_qc_summary(varus_stats_path)

    # Gene set statistics section
    gene_stats_path = None
    for search_dir in [os.path.join(outdir, "quality_control"), outdir]:
        p = os.path.join(search_dir, "gene_set_statistics.txt")
        if os.path.exists(p):
            gene_stats_path = p
            break

    if gene_stats_path:
        html += "\n<h2>Gene Set Statistics</h2>\n"
        html += read_qc_summary(gene_stats_path)

    # Gene set plots
    gene_plots = [
        ("isoform_and_exon_structure", "Isoform distribution (left) and single-exon vs multi-exon transcript structure (right)."),
        ("transcript_lengths", "CDS length excluding introns (left) and genomic span including introns (right)."),
        ("introns_per_gene", "Number of introns per gene."),
        ("evidence_support", "Evidence support from transcript (RNA-Seq or IsoSeq) and protein hints."),
    ]
    for plot_name, caption in gene_plots:
        for search_dir in [os.path.join(outdir, "quality_control"), outdir]:
            path = os.path.join(search_dir, f"{plot_name}.png")
            if os.path.exists(path):
                img_html = embed_image(path, download_name=f"{plot_name}.png")
                if img_html:
                    html += f'<div class="figure">\n{img_html}\n'
                    html += f'<div class="caption">{caption}</div>\n</div>\n'
                break

    # Training summary plot
    if "training_summary" in images and images["training_summary"]:
        # Detect dual mode from directory layout
        is_dual_etp = (os.path.isdir(os.path.join(outdir, "GeneMark-ETP")) and
                       os.path.isdir(os.path.join(outdir, "GeneMark-ETP-isoseq"))) or \
                      (os.path.isdir(os.path.join(outdir, "..", "GeneMark-ETP")) and
                       os.path.isdir(os.path.join(outdir, "..", "GeneMark-ETP-isoseq")))
        caption = (
            "Training gene counts at each pipeline stage (left) and AUGUSTUS parameter "
            "accuracy before and after optimization (right)."
        )
        if is_dual_etp:
            caption += (
                " In dual mode, GeneMark-ETP was run twice in parallel (short-read and IsoSeq); "
                "the training genes from both runs were merged and deduplicated before AUGUSTUS training."
            )
        html += f"""
<h2>Training Summary</h2>
<div class="figure">
{images["training_summary"]}
<div class="caption">{caption}</div>
</div>
"""

    # Completeness plot (BUSCO + compleasm combined)
    completeness_img = None
    for search_dir in [os.path.join(outdir, "quality_control"), outdir]:
        p = os.path.join(search_dir, "completeness.png")
        if os.path.exists(p):
            completeness_img = embed_image(p, download_name="completeness.png")
            break

    if completeness_img:
        html += """
<h2>Completeness Assessment</h2>
<div class="figure">
"""
        html += completeness_img
        html += """
<div class="caption">BUSCO and compleasm completeness assessment of the genome assembly and
predicted proteome. S=Single-copy, D=Duplicated, F=Fragmented, M=Missing.</div>
</div>
"""

    # FANTASIA functional annotation (optional, GPU-only)
    fantasia_summary_text = images.get("fantasia_summary_text") or ""
    fantasia_plot_html    = images.get("fantasia_plot") or ""
    if fantasia_summary_text or fantasia_plot_html:
        html += """
<h2>Functional Annotation (FANTASIA-Lite)</h2>
<p>Predicted proteins were functionally annotated with
<a href="https://doi.org/10.1038/s42003-025-08651-2" target="_blank">FANTASIA-Lite</a>,
which assigns Gene Ontology terms via ProtT5 protein language model embeddings
against a bundled lookup of pre-computed reference embeddings. This step is
optional and runs only when <code>fantasia.enable=1</code>.</p>
"""
        if fantasia_summary_text:
            html += "<pre>" + fantasia_summary_text + "</pre>\n"
        if fantasia_plot_html:
            html += '<div class="figure">\n'
            html += fantasia_plot_html
            html += '\n<div class="caption">FANTASIA-Lite functional category distribution. Each predicted protein is assigned to every curated functional category its GO terms touch (energy / photosynthesis, carbohydrate metabolism, transport, signal transduction, and so on); pie slices show the fraction of all category memberships in each bucket. Proteins matching no curated category fall into "Other / unclassified". The raw GO term assignments per namespace are still listed in the FANTASIA summary above.</div>\n</div>\n'

    # QC summaries
    if any(qc_data.values()):
        html += "\n<h2>Quality Control Summaries</h2>\n"
        for name, content in qc_data.items():
            if content:
                html += f"<h3>{name}</h3>\n{content}\n"

    # Runtime benchmarks
    if benchmarks:
        total_wall = sum(b["wallclock_s"] for b in benchmarks)
        total_cpu = sum(b["cpu_s"] for b in benchmarks)

        html += f"""
<h2>Runtime</h2>
<p>Total wall clock time: <strong>{format_time(total_wall)}</strong>.
Total CPU time: <strong>{format_time(total_cpu)}</strong>.</p>
"""

        # Runtime plot
        if "runtime_plot" in images and images["runtime_plot"]:
            html += '<div class="figure">\n'
            html += images["runtime_plot"]
            html += '\n<div class="caption">Resource consumption for the top 20 longest-running pipeline steps. Left: wall clock and CPU time. Right: peak RAM usage.</div>\n</div>\n'

        # Runtime table (collapsible)
        html += """
<details>
<summary>Full runtime table</summary>
<table>
<tr><th>Rule</th><th>Wall clock</th><th>CPU time</th><th>Max RSS (MB)</th></tr>
"""
        for b in benchmarks:
            if b["wallclock_s"] > 0.5:
                html += (
                    f'<tr><td>{b["rule"]}</td>'
                    f'<td>{format_time(b["wallclock_s"])}</td>'
                    f'<td>{format_time(b["cpu_s"])}</td>'
                    f'<td>{b["max_rss_mb"]:.0f}</td></tr>\n'
                )
        html += "</table>\n</details>\n"

    # Software Versions
    if software_versions:
        html += "\n<h2>Software Versions</h2>\n"
        html += "<p>The following software versions were used in this pipeline run:</p>\n"
        html += "<table>\n<tr><th>Software</th><th>Version</th></tr>\n"
        for tool, version in software_versions:
            html += f"<tr><td>{tool}</td><td>{version}</td></tr>\n"
        html += "</table>\n"

    # Citations
    if citations:
        html += """
<h2>References</h2>
<p>If you publish these BRAKER4 results, please cite the following publications:</p>
<ol class="citation-list">
"""
        for cite in citations:
            # Make DOI references clickable
            cite_html = re.sub(
                r'doi:(10\.\S+)',
                r'doi:<a href="https://doi.org/\1" target="_blank">\1</a>',
                cite
            )
            html += f"<li>{cite_html}</li>\n"
        html += "</ol>\n"

    html += f"""
<footer>
Generated by BRAKER4 v{pipeline_version} on {now}.
</footer>
</body>
</html>
"""
    return html


def main():
    parser = argparse.ArgumentParser(description="Generate BRAKER4 HTML report")
    parser.add_argument("-d", "--workdir", required=True,
                        help="Sample working directory (output/{sample}/)")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory for report files")
    parser.add_argument("-s", "--sample", default=None,
                        help="Sample name (auto-detected from workdir if omitted)")
    parser.add_argument("--mode", default=None,
                        help="Snakemake mode key (es/et/ep/etp/isoseq/dual); "
                             "overrides filesystem-based detect_mode()")
    args = parser.parse_args()

    d = args.workdir
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    sample_name = args.sample or os.path.basename(d.rstrip("/"))

    _MODE_LABELS = {
        "es":     "ES (ab initio)",
        "et":     "ET/BRAKER1 (RNA-Seq only)",
        "ep":     "EP/BRAKER2 (proteins only)",
        "etp":    "ETP/BRAKER3 (RNA-Seq + proteins)",
        "isoseq": "IsoSeq/BRAKER3 (long-read + proteins)",
        "dual":   "Dual/BRAKER3 (short-read + IsoSeq + proteins)",
    }
    if args.mode and args.mode in _MODE_LABELS:
        mode = _MODE_LABELS[args.mode]
    else:
        mode = detect_mode(d)

    # Generate methods text from pipeline outputs (single cohesive paragraph)
    methods_text = generate_methods_text(d, mode)

    # Read citations accumulated by rules
    citations_text = read_file(os.path.join(d, "report_citations.txt"))
    bibtex_text = read_file(os.path.join(d, "report_citations.bib"))

    # Deduplicate
    citations = deduplicate_citations(citations_text)
    bibtex_dedup = deduplicate_bibtex(bibtex_text)

    # Embed images
    images = {}
    for name, filename in [
        ("training_summary", "training_summary.png"),
    ]:
        # Check both QC subdir and results dir
        for search_dir in [os.path.join(out, "quality_control"), out, d]:
            path = os.path.join(search_dir, filename)
            if os.path.exists(path):
                images[name] = embed_image(path, download_name=filename)
                break

    # QC summaries
    qc_data = OrderedDict()
    for name, filename in [
        ("BUSCO", "busco_summary.txt"),
        ("Compleasm (Proteome)", "compleasm_summary.txt"),
        ("OMArk", "omark_summary.txt"),
        ("gffcompare", "gffcompare.stats"),
    ]:
        for search_dir in [os.path.join(out, "quality_control"), out, d]:
            path = os.path.join(search_dir, filename)
            if os.path.exists(path):
                is_compleasm = "compleasm" in name.lower()
                qc_data[name] = read_qc_summary(path, reformat_compleasm=is_compleasm)
                break

    # FANTASIA-Lite functional annotation (optional). The collect_results rule
    # places these in <results>/quality_control/fantasia/, but the workdir copy
    # is also checked for in-place report runs.
    for search_dir in [
        os.path.join(out, "quality_control", "fantasia"),
        os.path.join(d, "fantasia"),
    ]:
        plot_path = os.path.join(search_dir, "fantasia_go_categories.png")
        summary_path = os.path.join(search_dir, "fantasia_summary.txt")
        if os.path.exists(plot_path) or os.path.exists(summary_path):
            if os.path.exists(plot_path):
                images["fantasia_plot"] = embed_image(
                    plot_path, download_name="fantasia_go_categories.png"
                )
            if os.path.exists(summary_path):
                images["fantasia_summary_text"] = read_file(summary_path)
            break

    # Compleasm genome assessment (separate from protein assessment above)
    for search_dir in [d, os.path.join(out, "quality_control"), out]:
        genome_summary = os.path.join(search_dir, "compleasm_genome_out", "summary.txt")
        if os.path.exists(genome_summary):
            qc_data["Compleasm (Genome)"] = read_qc_summary(genome_summary, reformat_compleasm=True)
            # Move genome entry before proteome entry
            if "Compleasm (Proteome)" in qc_data:
                items = list(qc_data.items())
                genome_idx = next(i for i, (k, _) in enumerate(items) if k == "Compleasm (Genome)")
                proteome_idx = next(i for i, (k, _) in enumerate(items) if k == "Compleasm (Proteome)")
                if genome_idx > proteome_idx:
                    genome_item = items.pop(genome_idx)
                    items.insert(proteome_idx, genome_item)
                    qc_data = OrderedDict(items)
            break

    # Collect benchmarks
    benchmarks = collect_benchmarks(d)

    # Generate runtime plot and embed it
    if benchmarks:
        qc_dir = os.path.join(out, "quality_control")
        os.makedirs(qc_dir, exist_ok=True)
        plot_path = generate_benchmark_plot(benchmarks, qc_dir)
        if plot_path:
            images["runtime_plot"] = embed_image(plot_path, download_name="runtime_plot.png")

    # Read software versions
    software_versions = read_software_versions(out, d)

    # Generate HTML
    html = generate_html(sample_name, mode, methods_text, citations, images, qc_data, benchmarks, outdir=out, software_versions=software_versions)
    html_path = os.path.join(out, "braker_report.html")
    with open(html_path, "w") as f:
        f.write(html)
    print(f"HTML report: {html_path}")

    # Write deduplicated BibTeX
    if bibtex_dedup:
        bib_path = os.path.join(out, "braker_citations.bib")
        with open(bib_path, "w") as f:
            f.write(bibtex_dedup + "\n")
        print(f"BibTeX file: {bib_path}")


if __name__ == "__main__":
    main()
