"""
Prepare protein hints for GeneMark-EP from ProtHint output.

Extracts protein-based intron hints (src=P) from prothint_augustus.gff
and converts them to GeneMark-compatible format. Also creates the
evidence GFF file for GeneMark-EP+ mode (src=M hints).

Mirrors braker.pl's get_genemark_hints() for EP mode and
create_evidence_gff() for EP mode.

Input:
    - prothint_hints.gff - ProtHint output with protein hints

Output:
    - genemark_hints_ep.gff - Protein intron hints for GeneMark-EP
    - genemark_evidence.gff - Evidence hints (src=M) for GeneMark-EP+

Container: teambraker/braker3:latest (contains join_mult_hints.pl)
"""

def _get_prepare_ep_inputs(wildcards):
    """ProtHint hints are always required; compleasm hints are added as a
    second src=M source when use_compleasm_hints is on (the default).
    braker.pl extracts src=M hints from the merged hintsfile (which
    contains compleasm + manual hints). BRAKER4's order has
    prepare_genemark_hints_ep running BEFORE merge_hints, so we read
    directly from compleasm_hints.gff instead.
    """
    inputs = {"prothint": f"output/{wildcards.sample}/prothint_hints.gff"}
    if config.get('use_compleasm_hints', True):
        inputs["compleasm_hints"] = f"output/{wildcards.sample}/compleasm_hints.gff"
    return inputs


rule prepare_genemark_hints_ep:
    input:
        unpack(_get_prepare_ep_inputs)
    output:
        gm_hints="output/{sample}/genemark_hints_ep.gff",
        evidence="output/{sample}/genemark_evidence.gff"
    log:
        "logs/{sample}/prepare_genemark_hints_ep/prepare.log"
    benchmark:
        "benchmarks/{sample}/prepare_genemark_hints_ep/prepare_genemark_hints_ep.txt"
    params:
        # Resolved at scheduling time. Empty string when use_compleasm_hints
        # is off, so the shell loop just iterates over an empty placeholder.
        compleasm_hints_path = lambda w: (
            f"output/{w.sample}/compleasm_hints.gff"
            if config.get('use_compleasm_hints', True) else ""
        )
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "Extracting protein hints for GeneMark-EP..." > {log}

        # Extract src=P hints and convert to GeneMark format
        # - Convert feature types: Intron->intron, start->start_codon, stop->stop_codon
        # - Set score to multiplicity value
        # Using perl instead of awk for portability (match with capture groups)
        perl -F'\t' -lane '
            next unless /src=P/;
            $F[2] =~ s/^Intron$/intron/;
            $F[2] =~ s/^start$/start_codon/;
            $F[2] =~ s/^stop$/stop_codon/;
            if ($F[8] =~ /mult=([^;]+)/) {{
                $F[5] = $1;
            }} else {{
                $F[5] = 1;
            }}
            print join("\t", @F);
        ' {input.prothint} > {output.gm_hints}.unsorted

        if [ -s {output.gm_hints}.unsorted ]; then
            # Sort and join multiplicities
            cat {output.gm_hints}.unsorted | \
                sort -n -k 4,4 | sort -s -n -k 5,5 | \
                sort -s -k 3,3 | sort -s -k 1,1 | \
                join_mult_hints.pl > {output.gm_hints}
        else
            touch {output.gm_hints}
            echo "WARNING: No src=P hints found in ProtHint output" >> {log}
        fi

        rm -f {output.gm_hints}.unsorted

        n_hints=$(wc -l < {output.gm_hints})
        echo "Extracted $n_hints protein hints for GeneMark-EP" >> {log}

        # Create evidence GFF (src=M hints for GeneMark-EP+ mode).
        # src=M hints come from compleasm and/or user-supplied manual hints.
        # ProtHint output occasionally also contains src=M lines (start/stop
        # codon anchors), so we scan both sources and concatenate.
        # braker.pl extracts these from the merged hintsfile.gff after all hint
        # sources are combined; BRAKER4 reads directly from each source because
        # prepare_genemark_hints_ep runs before merge_hints in our DAG.
        > {output.evidence}
        for src_file in {input.prothint} "{params.compleasm_hints_path}"; do
            [ -n "$src_file" ] && [ -f "$src_file" ] || continue
            perl -F'\t' -lane '
                next unless /src=M/;
                $F[2] =~ s/^Intron$/intron/;
                $F[2] =~ s/^start$/start_codon/;
                $F[2] =~ s/^stop$/stop_codon/;
                print join("\t", @F);
            ' "$src_file" >> {output.evidence}
        done

        n_evidence=$(wc -l < {output.evidence})
        echo "Extracted $n_evidence evidence hints (src=M from prothint + compleasm if available)" >> {log}
        """
