# Graph Report - /home/katharina/git/BRAKER4  (2026-04-18)

## Corpus Check
- 17 files · ~1,098,877 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 389 nodes · 581 edges · 28 communities detected
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 37 edges (avg confidence: 0.82)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]

## God Nodes (most connected - your core abstractions)
1. `BRAKER4 Test Scenarios README` - 19 edges
2. `BRAKER4 pipeline overview diagram` - 19 edges
3. `main()` - 16 edges
4. `README.md (BRAKER4)` - 15 edges
5. `main()` - 12 edges
6. `main()` - 12 edges
7. `BRAKER4 pipeline overview diagram (SVG, Graphviz-generated)` - 12 edges
8. `test_data/genome.fa (A. thaliana chr5 fragment, 8 contigs)` - 10 edges
9. `Local Test Scenarios (scenarios 01-07, A. thaliana fragment)` - 9 edges
10. `Pipeline GeneMark node: ES/ET/EP+/ETP modes, auto mode selection` - 9 edges

## Surprising Connections (you probably didn't know these)
- `braker.gtf output file` --semantically_similar_to--> `Pipeline output: braker.gtf, braker.gff3, braker.aa, braker.codingseq`  [INFERRED] [semantically similar]
  test_scenarios/README.md → img/pipeline_overview.pdf
- `quality_control/ output directory (BUSCO, compleasm, OMArk)` --semantically_similar_to--> `Quality control step (BUSCO, compleasm, OMArk optional, gffcompare optional, gene-support)`  [INFERRED] [semantically similar]
  test_scenarios/README.md → img/pipeline_overview.pdf
- `ES mode (ab initio, no evidence)` --semantically_similar_to--> `GeneMark step (ES/ET/EP+/ETP, auto-mode selection)`  [INFERRED] [semantically similar]
  test_scenarios/README.md → img/pipeline_overview.pdf
- `[fantasia] section in HPC config.ini (FANTASIA-Lite settings)` --semantically_similar_to--> `Functional annotation step - FANTASIA-Lite (ProtT5 embeddings to GO terms, GPU only) - optional`  [INFERRED] [semantically similar]
  test_scenarios/config.ini → img/pipeline_overview.pdf
- `local scenario_08_ep_redmask samples.csv (EP mode, unmasked genome)` --references--> `Repeat masking step (RepeatModeler2 + RepeatMasker + TRF)`  [INFERRED]
  test_scenarios_local/scenario_08_ep_redmask/samples.csv → img/pipeline_overview.pdf

## Hyperedges (group relationships)
- **GeneMark trains AUGUSTUS, both merged by TSEBRA: core gene prediction flow** — concept_genemark, concept_augustus, concept_tsebra [EXTRACTED 1.00]
- **Mode auto-detection from samples.csv evidence columns (ET/EP/ETP/ES/IsoSeq/Dual)** — samples_csv, concept_et_mode, concept_ep_mode, concept_etp_mode, concept_es_mode, concept_isoseq_mode, concept_dual_mode [EXTRACTED 1.00]
- **IsoSeq FASTQ aligned with minimap2+minisplice inside dedicated container** — testdata_isoseq_fastq, concept_minimap2, concept_minisplice, dockerfile_minisplice [EXTRACTED 1.00]
- **Core BRAKER4 gene prediction flow: GeneMark trains AUGUSTUS which merges via TSEBRA** — pipeline_step_genemark, pipeline_step_augustus, pipeline_step_tsebra [EXTRACTED 1.00]
- **Test configuration system: compute_profile.sh, config.ini, and scenario_overrides.sh jointly define all test execution parameters** — readme_compute_profile_sh, readme_config_ini, readme_scenario_overrides_sh [EXTRACTED 1.00]
- **HPC scenarios 08-11 collectively cover all four GeneMark modes (ES, EP, ET, ETP) on O. tauri genome** — scenario_08_es_samples, scenario_09_ep_masking_samples, scenario_10_et_varus_samples, scenario_11_etp_sra_samples [EXTRACTED 1.00]
- **BRAKER4 brand identity: three logo variants (hatched, hatching, primary) all depict hammer + egg motif representing the tool name and theme** — img_hatched_logo, img_hatching_logo, img_logo [INFERRED 0.90]
- **Core gene prediction pipeline: evidence hints feed GeneMark and AUGUSTUS, which are merged by TSEBRA — the central required workflow in all scenarios** — pipeline_node_evidence_hints, pipeline_node_genemark, pipeline_node_augustus, pipeline_node_tsebra [EXTRACTED 1.00]
- **All scenario rulegraphs instantiate the same BRAKER4 Snakemake workflow with different input evidence types (ES, EP+, ET, ETP, IsoSeq, dual, multi-mode)** — scenario_08_es_rulegraph, scenario_09_ep_masking_rulegraph, scenario_10_et_varus_rulegraph, scenario_11_etp_sra_rulegraph, scenario_12_multi_mode_rulegraph, scenario_local_01_es_rulegraph, scenario_local_02_ep_rulegraph, scenario_local_03_et_fastq_rulegraph, scenario_local_04_etp_bam_rulegraph, scenario_local_05_isoseq_bam_rulegraph, scenario_local_06_isoseq_fastq_rulegraph, scenario_local_07_dual_rulegraph, scenario_local_08_ep_redmask_rulegraph [INFERRED 0.88]

## Communities

### Community 0 - "Community 0"
Cohesion: 0.07
Nodes (49): BUSCO lineage eukaryota_odb12, local scenario_01_es samples.csv (ES mode, A. thaliana), local scenario_02_ep samples.csv (EP mode, A. thaliana), local scenario_03_et_fastq samples.csv (ET mode, FASTQ, A. thaliana), local scenario_04_etp_bam samples.csv (ETP mode, BAM+proteins), local scenario_05_isoseq_bam samples.csv (IsoSeq BAM multi-lib + proteins), local scenario_06_isoseq_fastq samples.csv (IsoSeq FASTQ + proteins), local scenario_07_dual samples.csv (Dual mode: BAM + IsoSeq + proteins) (+41 more)

### Community 1 - "Community 1"
Cohesion: 0.08
Nodes (38): add_intron_features(), build_tree(), check_overlap_compatibility(), compute_utr_features(), construct_gene_tree(), construct_transcript_tree(), create_introns_hash(), extract_introns_from_dict() (+30 more)

### Community 2 - "Community 2"
Cohesion: 0.08
Nodes (38): AUGUSTUS (training + hint-guided prediction), best_by_compleasm BUSCO-rescue in TSEBRA merge, teambraker/braker3:v3.0.10 container image, compleasm CDSpart hints for AUGUSTUS, Dual mode (IsoSeq + short-read RNA-Seq + proteins), Environment variable overrides (BRAKER4_<KEY>), EP mode (proteins only / BRAKER2 equivalent), ES mode (genome only / ab initio) (+30 more)

### Community 3 - "Community 3"
Cohesion: 0.1
Nodes (31): collect_benchmarks(), deduplicate_bibtex(), deduplicate_citations(), detect_mode(), embed_image(), format_bbc_decisions(), format_compleasm_as_busco(), format_time() (+23 more)

### Community 4 - "Community 4"
Cohesion: 0.13
Nodes (26): check_binary(), check_dir(), check_file(), determine_mode(), find_genemark_gtf(), find_input_files(), load_length_cutoff(), load_score_cutoff() (+18 more)

### Community 5 - "Community 5"
Cohesion: 0.14
Nodes (24): BRAKER4 pipeline overview diagram (PNG raster), BRAKER4 pipeline overview diagram (SVG, Graphviz-generated), Pipeline AUGUSTUS node: training + prediction (etraining, optimize, hint-guided prediction), Pipeline evidence hints node: RNA-Seq, Proteins (ProtHint), IsoSeq, BUSCO, Pipeline functional annotation node (optional, GPU): FANTASIA-Lite, ProtT5 embeddings, GO terms, GFF3 Ontology_term decoration, Pipeline GeneMark node: ES/ET/EP+/ETP modes, auto mode selection, Pipeline INPUTS cluster: Short-read RNA-Seq, IsoSeq long reads, Genome FASTA, Protein FASTA, Pipeline ncRNA annotation node (optional): barrnap, tRNAscan-SE, Infernal+Rfam, FEELnc (+16 more)

### Community 6 - "Community 6"
Cohesion: 0.21
Nodes (20): Pipeline input: Genome FASTA, Pipeline input: IsoSeq long reads, Pipeline input: Protein FASTA, Pipeline input: Short-read RNA-Seq, Pipeline output: braker.gtf, braker.gff3, braker.aa, braker.codingseq, Pipeline output: braker_utr.gtf, Pipeline output: braker.go.gff3, fantasia_go_terms.tsv, Pipeline output: HTML report + citations.bib + software_versions.tsv (+12 more)

### Community 7 - "Community 7"
Cohesion: 0.16
Nodes (17): analyze_hint_support(), check_cds_overlap(), check_intron_support(), load_hints(), load_transcripts_and_introns(), main(), parse_gff_line(), parse_gtf_line() (+9 more)

### Community 8 - "Community 8"
Cohesion: 0.19
Nodes (15): get_cds_features(), get_sequence(), main(), normalize_transcript(), parse_attrs(), parse_gtf(), Update gene and transcript boundaries to match their features., Write normalized GTF, preserving gene/transcript/feature order. (+7 more)

### Community 9 - "Community 9"
Cohesion: 0.14
Nodes (10): check_tool_in_given_path(), create_log_file_name(), create_random_string(), create_tmp_dir(), find_tool(), Funtion that creates a random string added to the logfile name         and tmp d, Function that creates a log file with a random name, Function that creates a directory for temporary files with a random name (+2 more)

### Community 10 - "Community 10"
Cohesion: 0.27
Nodes (9): Enum, addToDict(), addToGc(), initCounts(), main(), parse(), parseCmd(), printStatistics() (+1 more)

### Community 11 - "Community 11"
Cohesion: 0.24
Nodes (11): compute_statistics(), generate_plots(), main(), parse_gtf(), parse_support_tsv(), Parse gene_support.tsv for evidence support visualization., Generate all publication-quality plots., Parse GTF file into gene/transcript/exon structure.      Returns:         genes: (+3 more)

### Community 12 - "Community 12"
Cohesion: 0.26
Nodes (11): classify_proteins(), _compile_categories(), main(), parse_results(), Map each protein to all categories its GO terms touch.      Returns a Counter of, Stream results.csv: collect counters AND per-row tuples for the TSV., Flat per-(transcript, GO term) TSV with human-readable GO names., Pie chart of broad functional categories.      Each protein is counted once per (+3 more)

### Community 13 - "Community 13"
Cohesion: 0.29
Nodes (9): extract_gene_blocks(), main(), parse_cds_intron_chains(), parse_intron_hints(), Rename gene and transcript IDs in GTF lines to avoid collisions.      Prefixes a, Parse hintsfile.gff and return a set of (chrom, start, end, strand) for intron h, Parse a GTF and return per-transcript CDS intron chains.      Returns:         c, Extract full gene blocks (all GTF lines) for given transcript IDs.      Returns (+1 more)

### Community 14 - "Community 14"
Cohesion: 0.29
Nodes (9): _add_ontology_term(), decorate(), load_go_assignments(), _lookup_key(), main(), Two-pass: build transcript->gene map, then rewrite the GFF3.      Both `mRNA` an, Return the FANTASIA lookup key for a transcript-equivalent row.      Prefers `tr, Map transcript_id -> set of GO IDs above the score threshold. (+1 more)

### Community 15 - "Community 15"
Cohesion: 0.29
Nodes (9): check_exon_support(), check_support(), main(), parse_gtf_transcripts(), parse_hints(), Check how many intervals are supported by hints.      For introns: exact match (, Check how many exons have any overlapping exon/CDS hints.      Uses overlap (not, Parse hints file into interval trees by (chrom, feature_type, source_class). (+1 more)

### Community 16 - "Community 16"
Cohesion: 0.36
Nodes (7): generate_plot(), main(), parse_busco_summary(), parse_compleasm_summary(), Parse BUSCO summary for genome and proteome scores., Parse compleasm summary.txt., Generate horizontal stacked bar chart.

### Community 17 - "Community 17"
Cohesion: 0.32
Nodes (7): count_loci_in_gb(), main(), parse_accuracy_file(), parse_gene_count(), Parse accuracy file, return dict with nu_sen, nu_sp, ex_sen, ex_sp, gene_sen, ge, Extract a gene count from a file. Returns int or None., Count LOCUS entries in GenBank file.

### Community 18 - "Community 18"
Cohesion: 0.29
Nodes (7): config.ini HPC shared biology config, [fantasia] section in HPC config.ini (FANTASIA-Lite settings), [PARAMS] section in HPC config.ini, config.ini local shared biology config, [PARAMS] section in local config.ini, HuggingFace cache directory for FANTASIA offline mode, FANTASIA-Lite Singularity container (fantasia_lite.sif)

### Community 19 - "Community 19"
Cohesion: 0.6
Nodes (5): extract_tx_ids_from_tsv(), main(), miniprot_to_hints(), read_and_filter_gff(), run_simple_process()

### Community 20 - "Community 20"
Cohesion: 0.6
Nodes (4): main(), parse_tblout(), Parse Infernal --fmt 2 tblout output.      Fields (fmt 2):     0: idx, 1: target, write_gff3()

### Community 21 - "Community 21"
Cohesion: 1.0
Nodes (3): BRAKER4 hatched logo (chick with hammer, post-hatching), BRAKER4 hatching logo (hammer bursting from egg), BRAKER4 primary logo (hammer smashing decorated egg)

### Community 22 - "Community 22"
Cohesion: 1.0
Nodes (2): Translation table validation (codes 1, 6, 29 only), Rationale: restrict translation_table to {1,6,29} because gmes_petap.pl rejects others

### Community 23 - "Community 23"
Cohesion: 1.0
Nodes (2): Local scenario_05 rulegraph (IsoSeq BAM input): check_nooverlap_bam, merge_isoseq_bams, run_genemark_etp, compute_flanking_region, copy_augustus_config, create_new_species, convert_to_genbank, filter_genes_straining, redundancy_removal, downsample_training_genes, split_training_set, train_augustus, optimize_augustus, run_augustus_hints, merge_hints, get_anno_hints, fix_in_frame_stop_codons, run_tsebra, best_by_compleasm, both_agree_rescue, filter_internal_stop_codons, normalize_cds, add_utrs, gene_support_summary, fix_gtf, extract_final_sequences, get_longest_isoform, busco_proteins, busco_summary, access_completeness, convert_gtf_to_gff3, collect_results, all, Local scenario_06 rulegraph (IsoSeq FASTQ input): prepare_masked_genome, minimap2_isoseq_align, sort_isoseq_bam, run_genemark_etp, compute_flanking_region, copy_augustus_config, create_new_species, convert_to_genbank, filter_genes_straining, redundancy_removal, downsample_training_genes, split_training_set, train_augustus, optimize_augustus, run_augustus_hints, merge_hints, get_anno_hints, fix_in_frame_stop_codons, run_tsebra, best_by_compleasm, both_agree_rescue, filter_internal_stop_codons, normalize_cds, add_utrs, gene_support_summary, fix_gtf, extract_final_sequences, get_longest_isoform, busco_proteins, busco_summary, access_completeness, convert_gtf_to_gff3, collect_results, all

### Community 24 - "Community 24"
Cohesion: 1.0
Nodes (0): 

### Community 25 - "Community 25"
Cohesion: 1.0
Nodes (1): .gitignore (BRAKER4)

### Community 26 - "Community 26"
Cohesion: 1.0
Nodes (1): MIT License (BRAKER4)

### Community 27 - "Community 27"
Cohesion: 1.0
Nodes (1): samples.csv + config.ini input specification

## Knowledge Gaps
- **124 isolated node(s):** `Parse BUSCO summary for genome and proteome scores.`, `Parse compleasm summary.txt.`, `Generate horizontal stacked bar chart.`, `Parse hintsfile.gff and return a set of (chrom, start, end, strand) for intron h`, `Parse a GTF and return per-transcript CDS intron chains.      Returns:         c` (+119 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 22`** (2 nodes): `Translation table validation (codes 1, 6, 29 only)`, `Rationale: restrict translation_table to {1,6,29} because gmes_petap.pl rejects others`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 23`** (2 nodes): `Local scenario_05 rulegraph (IsoSeq BAM input): check_nooverlap_bam, merge_isoseq_bams, run_genemark_etp, compute_flanking_region, copy_augustus_config, create_new_species, convert_to_genbank, filter_genes_straining, redundancy_removal, downsample_training_genes, split_training_set, train_augustus, optimize_augustus, run_augustus_hints, merge_hints, get_anno_hints, fix_in_frame_stop_codons, run_tsebra, best_by_compleasm, both_agree_rescue, filter_internal_stop_codons, normalize_cds, add_utrs, gene_support_summary, fix_gtf, extract_final_sequences, get_longest_isoform, busco_proteins, busco_summary, access_completeness, convert_gtf_to_gff3, collect_results, all`, `Local scenario_06 rulegraph (IsoSeq FASTQ input): prepare_masked_genome, minimap2_isoseq_align, sort_isoseq_bam, run_genemark_etp, compute_flanking_region, copy_augustus_config, create_new_species, convert_to_genbank, filter_genes_straining, redundancy_removal, downsample_training_genes, split_training_set, train_augustus, optimize_augustus, run_augustus_hints, merge_hints, get_anno_hints, fix_in_frame_stop_codons, run_tsebra, best_by_compleasm, both_agree_rescue, filter_internal_stop_codons, normalize_cds, add_utrs, gene_support_summary, fix_gtf, extract_final_sequences, get_longest_isoform, busco_proteins, busco_summary, access_completeness, convert_gtf_to_gff3, collect_results, all`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 24`** (1 nodes): `getAnnoFastaFromJoingenes.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 25`** (1 nodes): `.gitignore (BRAKER4)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 26`** (1 nodes): `MIT License (BRAKER4)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 27`** (1 nodes): `samples.csv + config.ini input specification`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ES mode (ab initio, no evidence)` connect `Community 0` to `Community 6`?**
  _High betweenness centrality (0.007) - this node is a cross-community bridge._
- **What connects `Parse BUSCO summary for genome and proteome scores.`, `Parse compleasm summary.txt.`, `Generate horizontal stacked bar chart.` to the rest of the system?**
  _124 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.07 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.08 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.08 - nodes in this community are weakly interconnected._
- **Should `Community 3` be split into smaller, more focused modules?**
  _Cohesion score 0.1 - nodes in this community are weakly interconnected._
- **Should `Community 4` be split into smaller, more focused modules?**
  _Cohesion score 0.13 - nodes in this community are weakly interconnected._