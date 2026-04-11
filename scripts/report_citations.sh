#!/bin/bash
# Central citation database for BRAKER4 report generation.
# Rules source this file and call: cite <key> <report_dir>
#
# Usage in a rule's shell block:
#   source {script_dir}/report_citations.sh
#   cite braker4 output/{wildcards.sample}
#   cite genemark_et output/{wildcards.sample}

cite() {
    local key="$1"
    local report_dir="$2"
    local txt_file="$report_dir/report_citations.txt"
    local bib_file="$report_dir/report_citations.bib"

    mkdir -p "$report_dir"

    (
    if command -v flock &>/dev/null; then flock 9; fi

    case "$key" in
        braker1)
            echo "Hoff, K. J., Lange, S., Lomsadze, A., Borodovsky, M., & Stanke, M. (2016). BRAKER1: unsupervised RNA-Seq-based genome annotation with GeneMark-ET and AUGUSTUS. Bioinformatics, 32(5), 767-769. doi:10.1093/bioinformatics/btv661" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{braker1,
  author  = {Hoff, Katharina J. and Lange, Simone and Lomsadze, Alexandre and Borodovsky, Mark and Stanke, Mario},
  title   = {{BRAKER1}: unsupervised {RNA-Seq}-based genome annotation with {GeneMark-ET} and {AUGUSTUS}},
  journal = {Bioinformatics},
  volume  = {32},
  number  = {5},
  pages   = {767--769},
  year    = {2016},
  doi     = {10.1093/bioinformatics/btv661}
}
BIBEOF
            ;;
        braker2)
            echo "Bruna, T., Hoff, K. J., Lomsadze, A., Stanke, M., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR Genomics and Bioinformatics, 3(1), lqaa108. doi:10.1093/nargab/lqaa108" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{braker2,
  author  = {Bruna, Tomas and Hoff, Katharina J. and Lomsadze, Alexandre and Stanke, Mario and Borodovsky, Mark},
  title   = {{BRAKER2}: Automatic Eukaryotic Genome Annotation with {GeneMark-EP+} and {AUGUSTUS} Supported by a Protein Database},
  journal = {NAR Genomics and Bioinformatics},
  volume  = {3},
  number  = {1},
  pages   = {lqaa108},
  year    = {2021},
  doi     = {10.1093/nargab/lqaa108}
}
BIBEOF
            ;;
        braker_book)
            echo "Bruna, T., Gabriel, L., & Hoff, K. J. (2025). Navigating Eukaryotic Genome Annotation Pipelines: A Route Map to Using BRAKER, Galba, and TSEBRA. Methods in Molecular Biology, 2935, 67-107. doi:10.1007/978-1-0716-4583-3_4" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@incollection{braker_book,
  author    = {Br{\r u}na, Tom{\'a}{\v s} and Gabriel, Lars and Hoff, Katharina J.},
  title     = {Navigating Eukaryotic Genome Annotation Pipelines: A Route Map to Using {BRAKER}, {Galba}, and {TSEBRA}},
  booktitle = {Methods in Molecular Biology},
  volume    = {2935},
  pages     = {67--107},
  year      = {2025},
  doi       = {10.1007/978-1-0716-4583-3_4}
}
BIBEOF
            ;;
        braker3)
            echo "Gabriel, L., Bruna, T., Hoff, K. J., Ebel, M., Lomsadze, A., Borodovsky, M., & Stanke, M. (2024). BRAKER3: Fully automated genome annotation using RNA-Seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA. Genome Research, 34(5), 769-777. doi:10.1101/gr.278090.123" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{braker3,
  author  = {Gabriel, Lars and Br{\r u}na, Tom{\'a}{\v s} and Hoff, Katharina J. and Ebel, Matthis and Lomsadze, Alexandre and Borodovsky, Mark and Stanke, Mario},
  title   = {{BRAKER3}: Fully automated genome annotation using {RNA-Seq} and protein evidence with {GeneMark-ETP}, {AUGUSTUS}, and {TSEBRA}},
  journal = {Genome Research},
  volume  = {34},
  number  = {5},
  pages   = {769--777},
  year    = {2024},
  doi     = {10.1101/gr.278090.123}
}
BIBEOF
            ;;
        genemark_es)
            echo "Lomsadze, A., Ter-Hovhannisyan, V., Chernoff, Y. O., & Borodovsky, M. (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. Nucleic Acids Research, 33(20), 6494-6506. doi:10.1093/nar/gki937" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{genemark_es,
  author  = {Lomsadze, Alexandre and Ter-Hovhannisyan, Vardges and Chernoff, Yury O. and Borodovsky, Mark},
  title   = {Gene identification in novel eukaryotic genomes by self-training algorithm},
  journal = {Nucleic Acids Research},
  volume  = {33},
  number  = {20},
  pages   = {6494--6506},
  year    = {2005},
  doi     = {10.1093/nar/gki937}
}
BIBEOF
            ;;
        genemark_et)
            echo "Lomsadze, A., Burns, P. D., & Borodovsky, M. (2014). Integration of mapped RNA-Seq reads into automatic training of eukaryotic gene finding algorithm. Nucleic Acids Research, 42(15), e119. doi:10.1093/nar/gku557" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{genemark_et,
  author  = {Lomsadze, Alexandre and Burns, Paul D. and Borodovsky, Mark},
  title   = {Integration of mapped {RNA-Seq} reads into automatic training of eukaryotic gene finding algorithm},
  journal = {Nucleic Acids Research},
  volume  = {42},
  number  = {15},
  pages   = {e119},
  year    = {2014},
  doi     = {10.1093/nar/gku557}
}
BIBEOF
            ;;
        genemark_ep)
            echo "Bruna, T., Lomsadze, A., & Borodovsky, M. (2020). GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins. NAR Genomics and Bioinformatics, 2(2), lqaa026. doi:10.1093/nargab/lqaa026" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{genemark_ep,
  author  = {Bruna, Tomas and Lomsadze, Alexandre and Borodovsky, Mark},
  title   = {{GeneMark-EP+}: eukaryotic gene prediction with self-training in the space of genes and proteins},
  journal = {NAR Genomics and Bioinformatics},
  volume  = {2},
  number  = {2},
  pages   = {lqaa026},
  year    = {2020},
  doi     = {10.1093/nargab/lqaa026}
}
BIBEOF
            ;;
        genemark_etp)
            echo "Bruna, T., Lomsadze, A., & Borodovsky, M. (2024). GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes. Genome Research, 34(5), 757-768. doi:10.1101/gr.278373.123" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{genemark_etp,
  author  = {Bruna, Tomas and Lomsadze, Alexandre and Borodovsky, Mark},
  title   = {{GeneMark-ETP} significantly improves the accuracy of automatic annotation of large eukaryotic genomes},
  journal = {Genome Research},
  volume  = {34},
  number  = {5},
  pages   = {757--768},
  year    = {2024},
  doi     = {10.1101/gr.278373.123}
}
BIBEOF
            ;;
        genemarks_t)
            echo "Tang, S., Lomsadze, A., & Borodovsky, M. (2015). Identification of protein coding regions in RNA transcripts. Nucleic Acids Research, 43(12), e78. doi:10.1093/nar/gkv227" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{genemarks_t,
  author  = {Tang, Shiyuyun and Lomsadze, Alexandre and Borodovsky, Mark},
  title   = {Identification of protein coding regions in {RNA} transcripts},
  journal = {Nucleic Acids Research},
  volume  = {43},
  number  = {12},
  pages   = {e78},
  year    = {2015},
  doi     = {10.1093/nar/gkv227}
}
BIBEOF
            ;;
        augustus)
            echo "Stanke, M., Diekhans, M., Baertsch, R., & Haussler, D. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics, 24(5), 637-644. doi:10.1093/bioinformatics/btn013" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{augustus,
  author  = {Stanke, Mario and Diekhans, Mark and Baertsch, Robert and Haussler, David},
  title   = {Using native and syntenically mapped {cDNA} alignments to improve de novo gene finding},
  journal = {Bioinformatics},
  volume  = {24},
  number  = {5},
  pages   = {637--644},
  year    = {2008},
  doi     = {10.1093/bioinformatics/btn013}
}
BIBEOF
            ;;
        tsebra)
            echo "Gabriel, L., Hoff, K. J., Bruna, T., Borodovsky, M., & Stanke, M. (2021). TSEBRA: transcript selector for BRAKER. BMC Bioinformatics, 22, 566. doi:10.1186/s12859-021-04482-0" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{tsebra,
  author  = {Gabriel, Lars and Hoff, Katharina J. and Bruna, Tomas and Borodovsky, Mark and Stanke, Mario},
  title   = {{TSEBRA}: transcript selector for {BRAKER}},
  journal = {BMC Bioinformatics},
  volume  = {22},
  pages   = {566},
  year    = {2021},
  doi     = {10.1186/s12859-021-04482-0}
}
BIBEOF
            ;;
        diamond)
            echo "Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18, 366-368. doi:10.1038/s41592-021-01101-x" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{diamond,
  author  = {Buchfink, Benjamin and Reuter, Klaus and Drost, Hajk-Georg},
  title   = {Sensitive protein alignments at tree-of-life scale using {DIAMOND}},
  journal = {Nature Methods},
  volume  = {18},
  pages   = {366--368},
  year    = {2021},
  doi     = {10.1038/s41592-021-01101-x}
}
BIBEOF
            ;;
        samtools)
            echo "Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi:10.1093/gigascience/giab008" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{samtools,
  author  = {Danecek, Petr and others},
  title   = {Twelve years of {SAMtools} and {BCFtools}},
  journal = {GigaScience},
  volume  = {10},
  number  = {2},
  pages   = {giab008},
  year    = {2021},
  doi     = {10.1093/gigascience/giab008}
}
BIBEOF
            ;;
        hisat2)
            echo "Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907-915. doi:10.1038/s41587-019-0201-4" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{hisat2,
  author  = {Kim, Daehwan and Paggi, Joseph M. and Park, Chanhee and Bennett, Christopher and Salzberg, Steven L.},
  title   = {Graph-based genome alignment and genotyping with {HISAT2} and {HISAT}-genotype},
  journal = {Nature Biotechnology},
  volume  = {37},
  number  = {8},
  pages   = {907--915},
  year    = {2019},
  doi     = {10.1038/s41587-019-0201-4}
}
BIBEOF
            ;;
        minimap2)
            echo "Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100. doi:10.1093/bioinformatics/bty191" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{minimap2,
  author  = {Li, Heng},
  title   = {Minimap2: pairwise alignment for nucleotide sequences},
  journal = {Bioinformatics},
  volume  = {34},
  number  = {18},
  pages   = {3094--3100},
  year    = {2018},
  doi     = {10.1093/bioinformatics/bty191}
}
BIBEOF
            ;;
        stringtie)
            echo "Kovaka, S., Zimin, A. V., Pertea, G. M., Razaghi, R., Salzberg, S. L., & Pertea, M. (2019). Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biology, 20, 278. doi:10.1186/s13059-019-1910-1" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{stringtie,
  author  = {Kovaka, Sam and Zimin, Aleksey V. and Pertea, Geo M. and Razaghi, Ryan and Salzberg, Steven L. and Pertea, Mihaela},
  title   = {Transcriptome assembly from long-read {RNA-seq} alignments with {StringTie2}},
  journal = {Genome Biology},
  volume  = {20},
  pages   = {278},
  year    = {2019},
  doi     = {10.1186/s13059-019-1910-1}
}
BIBEOF
            ;;
        busco)
            echo "Manni, M., Berkeley, M. R., Seppey, M., Simao, F. A., & Zdobnov, E. M. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, 38(10), 4647-4654. doi:10.1093/molbev/msab199" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{busco,
  author  = {Manni, Mosè and Berkeley, Matthew R. and Seppey, Mathieu and Simão, Felipe A. and Zdobnov, Evgeny M.},
  title   = {{BUSCO} Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes},
  journal = {Molecular Biology and Evolution},
  volume  = {38},
  number  = {10},
  pages   = {4647--4654},
  year    = {2021},
  doi     = {10.1093/molbev/msab199}
}
BIBEOF
            ;;
        compleasm)
            echo "Huang, N., & Li, H. (2023). compleasm: a faster and more accurate reimplementation of BUSCO. Bioinformatics, 39(10), btad595. doi:10.1093/bioinformatics/btad595" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{compleasm,
  author  = {Huang, Neng and Li, Heng},
  title   = {compleasm: a faster and more accurate reimplementation of {BUSCO}},
  journal = {Bioinformatics},
  volume  = {39},
  number  = {10},
  pages   = {btad595},
  year    = {2023},
  doi     = {10.1093/bioinformatics/btad595}
}
BIBEOF
            ;;
        omark)
            echo "Nevers, Y., Warwick Vesztrocy, A., Rossier, V., et al. (2025). Quality assessment of gene repertoire annotations with OMArk. Nature Biotechnology, 43, 124-133. doi:10.1038/s41587-024-02147-w" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{omark,
  author  = {Nevers, Yannis and Warwick Vesztrocy, Alex and Rossier, Victor and others},
  title   = {Quality assessment of gene repertoire annotations with {OMArk}},
  journal = {Nature Biotechnology},
  volume  = {43},
  pages   = {124--133},
  year    = {2025},
  doi     = {10.1038/s41587-024-02147-w}
}
BIBEOF
            ;;
        repeatmodeler)
            echo "Flynn, J. M., et al. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. PNAS, 117(17), 9451-9457. doi:10.1073/pnas.1921046117" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{repeatmodeler,
  author  = {Flynn, Jullien M. and others},
  title   = {{RepeatModeler2} for automated genomic discovery of transposable element families},
  journal = {PNAS},
  volume  = {117},
  number  = {17},
  pages   = {9451--9457},
  year    = {2020},
  doi     = {10.1073/pnas.1921046117}
}
BIBEOF
            ;;
        agat)
            echo "Dainat, J. (2024). AGAT: Another Gff Analysis Toolkit. doi:10.5281/zenodo.3552717" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@software{agat,
  author  = {Dainat, Jacques},
  title   = {{AGAT}: Another Gff Analysis Toolkit},
  year    = {2024},
  doi     = {10.5281/zenodo.3552717}
}
BIBEOF
            ;;
        barrnap)
            echo "Seemann, T. (2018). barrnap: BAsic Rapid Ribosomal RNA Predictor. https://github.com/tseemann/barrnap" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@software{barrnap,
  author  = {Seemann, Torsten},
  title   = {barrnap: {BAsic} Rapid Ribosomal {RNA} Predictor},
  year    = {2018},
  url     = {https://github.com/tseemann/barrnap}
}
BIBEOF
            ;;
        varus)
            echo "Stanke, M., Bruhn, W., Becker, F., & Hoff, K. J. (2019). VARUS: sampling complementary RNA reads from the Sequence Read Archive. BMC Bioinformatics, 20, 558. doi:10.1186/s12859-019-3182-x" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{varus,
  author  = {Stanke, Mario and Bruhn, Willy and Becker, Felix and Hoff, Katharina J.},
  title   = {{VARUS}: sampling complementary {RNA} reads from the {Sequence Read Archive}},
  journal = {BMC Bioinformatics},
  volume  = {20},
  pages   = {558},
  year    = {2019},
  doi     = {10.1186/s12859-019-3182-x}
}
BIBEOF
            ;;
        prothint)
            echo "Bruna, T., Lomsadze, A., & Borodovsky, M. (2020). GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins. NAR Genomics and Bioinformatics, 2(2), lqaa026. doi:10.1093/nargab/lqaa026" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{prothint,
  author  = {Bruna, Tomas and Lomsadze, Alexandre and Borodovsky, Mark},
  title   = {{GeneMark-EP+}: eukaryotic gene prediction with self-training in the space of genes and proteins},
  journal = {NAR Genomics and Bioinformatics},
  volume  = {2},
  number  = {2},
  pages   = {lqaa026},
  year    = {2020},
  doi     = {10.1093/nargab/lqaa026}
}
BIBEOF
            ;;
        gffcompare)
            echo "Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. doi:10.12688/f1000research.23297.2" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{gffcompare,
  author  = {Pertea, Geo and Pertea, Mihaela},
  title   = {{GFF} Utilities: {GffRead} and {GffCompare}},
  journal = {F1000Research},
  volume  = {9},
  pages   = {304},
  year    = {2020},
  doi     = {10.12688/f1000research.23297.2}
}
BIBEOF
            ;;
        sratoolkit)
            echo "SRA Toolkit Development Team (2020). SRA Toolkit. https://github.com/ncbi/sra-tools" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@software{sratoolkit,
  author  = {{SRA Toolkit Development Team}},
  title   = {{SRA} Toolkit},
  year    = {2020},
  url     = {https://github.com/ncbi/sra-tools}
}
BIBEOF
            ;;
        repeatmasker)
            echo "Smit, A. F. A., Hubley, R., & Green, P. (2013). RepeatMasker Open-4.0. http://www.repeatmasker.org" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@software{repeatmasker,
  author  = {Smit, A. F. A. and Hubley, Robert and Green, Phil},
  title   = {{RepeatMasker} Open-4.0},
  year    = {2013},
  url     = {http://www.repeatmasker.org}
}
BIBEOF
            ;;
        trf)
            echo "Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. Nucleic Acids Research, 27(2), 573-580. doi:10.1093/nar/27.2.573" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{trf,
  author  = {Benson, Gary},
  title   = {Tandem repeats finder: a program to analyze {DNA} sequences},
  journal = {Nucleic Acids Research},
  volume  = {27},
  number  = {2},
  pages   = {573--580},
  year    = {1999},
  doi     = {10.1093/nar/27.2.573}
}
BIBEOF
            ;;
        spaln)
            echo "Iwata, H., & Gotoh, O. (2012). Benchmarking spliced alignment programs including Spaln2, an extended version of Spaln that incorporates additional species-specific features. Nucleic Acids Research, 40(20), e161. doi:10.1093/nar/gks708" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{spaln,
  author  = {Iwata, Hiroaki and Gotoh, Osamu},
  title   = {Benchmarking spliced alignment programs including {Spaln2}, an extended version of {Spaln} that incorporates additional species-specific features},
  journal = {Nucleic Acids Research},
  volume  = {40},
  number  = {20},
  pages   = {e161},
  year    = {2012},
  doi     = {10.1093/nar/gks708}
}
BIBEOF
            ;;
        miniprot)
            echo "Li, H. (2023). Protein-to-genome alignment with miniprot. Bioinformatics, 39(1), btad014. doi:10.1093/bioinformatics/btad014" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{miniprot,
  author  = {Li, Heng},
  title   = {Protein-to-genome alignment with miniprot},
  journal = {Bioinformatics},
  volume  = {39},
  number  = {1},
  pages   = {btad014},
  year    = {2023},
  doi     = {10.1093/bioinformatics/btad014}
}
BIBEOF
            ;;
        omamer)
            echo "Rossier, V., Warwick Vesztrocy, A., Robinson-Rechavi, M., & Dessimoz, C. (2021). OMAmer: tree-driven and alignment-free protein assignment to subfamilies outperforms closest sequence approaches. Bioinformatics, 37(18), 2866-2873. doi:10.1093/bioinformatics/btab219" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{omamer,
  author  = {Rossier, Victor and Warwick Vesztrocy, Alex and Robinson-Rechavi, Marc and Dessimoz, Christophe},
  title   = {{OMAmer}: tree-driven and alignment-free protein assignment to subfamilies outperforms closest sequence approaches},
  journal = {Bioinformatics},
  volume  = {37},
  number  = {18},
  pages   = {2866--2873},
  year    = {2021},
  doi     = {10.1093/bioinformatics/btab219}
}
BIBEOF
            ;;
        snakemake)
            echo "Molder, F., et al. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33. doi:10.12688/f1000research.29032.2" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{snakemake,
  author  = {Mölder, Felix and others},
  title   = {Sustainable data analysis with {Snakemake}},
  journal = {F1000Research},
  volume  = {10},
  pages   = {33},
  year    = {2021},
  doi     = {10.12688/f1000research.29032.2}
}
BIBEOF
            ;;
        braker4)
            echo "Krall, H., & Hoff, K. J. (2026). BRAKER4: a Snakemake pipeline for fully automated eukaryotic genome annotation. In preparation." >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@unpublished{braker4,
  author  = {Krall, Henning and Hoff, Katharina J.},
  title   = {{BRAKER4}: a {Snakemake} pipeline for fully automated eukaryotic genome annotation},
  year    = {2026},
  note    = {In preparation}
}
BIBEOF
            ;;
        trnascan)
            echo "Chan, P. P., & Lowe, T. M. (2019). tRNAscan-SE: searching for tRNA genes in genomic sequences. In Gene Prediction (pp. 1-14). Humana, New York, NY. doi:10.1007/978-1-4939-9173-0_1" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@incollection{trnascan,
  author    = {Chan, Patricia P. and Lowe, Todd M.},
  title     = {{tRNAscan-SE}: searching for {tRNA} genes in genomic sequences},
  booktitle = {Gene Prediction},
  pages     = {1--14},
  publisher = {Humana, New York, NY},
  year      = {2019},
  doi       = {10.1007/978-1-4939-9173-0_1}
}
BIBEOF
            ;;
        infernal)
            echo "Nawrocki, E. P., & Eddy, S. R. (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29(22), 2933-2935. doi:10.1093/bioinformatics/btt509" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{infernal,
  author  = {Nawrocki, Eric P. and Eddy, Sean R.},
  title   = {Infernal 1.1: 100-fold faster {RNA} homology searches},
  journal = {Bioinformatics},
  volume  = {29},
  number  = {22},
  pages   = {2933--2935},
  year    = {2013},
  doi     = {10.1093/bioinformatics/btt509}
}
BIBEOF
            ;;
        rfam)
            echo "Kalvari, I., Nawrocki, E. P., Ontiveros-Palacios, N., et al. (2021). Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Research, 49(D1), D192-D200. doi:10.1093/nar/gkaa1047" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{rfam,
  author  = {Kalvari, Ioanna and Nawrocki, Eric P. and Ontiveros-Palacios, Nancy and others},
  title   = {Rfam 14: expanded coverage of metagenomic, viral and {microRNA} families},
  journal = {Nucleic Acids Research},
  volume  = {49},
  number  = {D1},
  pages   = {D192--D200},
  year    = {2021},
  doi     = {10.1093/nar/gkaa1047}
}
BIBEOF
            ;;
        fantasia)
            echo "Martínez-Redondo, G. I., Perez-Canales, F. M., Carbonetto, B., Fernández, J. M., Barrios-Núñez, I., Vázquez-Valls, M., Cases, I., Rojas, A. M., & Fernández, R. (2025). FANTASIA leverages language models to decode the functional dark proteome across the animal tree of life. Communications Biology, 8, 1227. doi:10.1038/s42003-025-08651-2" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{fantasia,
  author  = {Mart{\'\i}nez-Redondo, Gemma I. and Perez-Canales, Francisco M. and Carbonetto, Bel{\'e}n and Fern{\'a}ndez, Jos{\'e} M. and Barrios-N{\'u}{\~n}ez, Israel and V{\'a}zquez-Valls, Mar{\c c}al and Cases, Ildefonso and Rojas, Ana M. and Fern{\'a}ndez, Rosa},
  title   = {{FANTASIA} leverages language models to decode the functional dark proteome across the animal tree of life},
  journal = {Communications Biology},
  volume  = {8},
  pages   = {1227},
  year    = {2025},
  doi     = {10.1038/s42003-025-08651-2}
}
BIBEOF
            ;;
        fantasia_methods)
            echo "Cases, I., Martínez-Redondo, G. I., Fernández, R., & Rojas, A. M. (2025). Functional Annotation of Proteomes Using Protein Language Models: A High-Throughput Implementation of the ProtTrans Model. Methods in Molecular Biology. doi:10.1007/978-1-0716-4623-6_8" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@incollection{fantasia_methods,
  author    = {Cases, Ildefonso and Mart{\'\i}nez-Redondo, Gemma I. and Fern{\'a}ndez, Rosa and Rojas, Ana M.},
  title     = {Functional Annotation of Proteomes Using Protein Language Models: A High-Throughput Implementation of the {ProtTrans} Model},
  booktitle = {Methods in Molecular Biology},
  year      = {2025},
  doi       = {10.1007/978-1-0716-4623-6_8}
}
BIBEOF
            ;;
        feelnc)
            echo "Wucher, V., Legeai, F., Hedan, B., et al. (2017). FEELnc: a tool for long non-coding RNA annotation and its application to the dog transcriptome. Nucleic Acids Research, 45(8), e57. doi:10.1093/nar/gkw1306" >> "$txt_file"
            cat >> "$bib_file" << 'BIBEOF'
@article{feelnc,
  author  = {Wucher, Valentin and Legeai, Fabrice and Hedan, Benoit and others},
  title   = {{FEELnc}: a tool for long non-coding {RNA} annotation and its application to the dog transcriptome},
  journal = {Nucleic Acids Research},
  volume  = {45},
  number  = {8},
  pages   = {e57},
  year    = {2017},
  doi     = {10.1093/nar/gkw1306}
}
BIBEOF
            ;;
        *)
            echo "WARNING: Unknown citation key: $key" >&2
            ;;
    esac

    ) 9>"$report_dir/.citations.lock"
}
