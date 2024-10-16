import pathlib

import pandas as pd

import snake_helpers

# read file overview: sample name, path to r1, path to r2, extra (ONT reads)
input_data = pd.read_csv(config["infile"], sep="\t")

ILM_ASSEMBLERS = {"skesa"}
HYBRID_ASSEMBLERS = {"sf_unicycler", "lf_polypolish"}
ALL_ASSEMBLERS = ILM_ASSEMBLERS.union(HYBRID_ASSEMBLERS)
SEROTYPERS = {"seqsero2", "sistr"}

rule all:
    input:
        all_serotyped = expand("typing/{serotyper}/{sample}_{assembler}_result.txt", # we need to move the result
        serotyper = SEROTYPERS, sample=input_data["sample"], assembler=ALL_ASSEMBLERS)

# we're letting shovill do the trimming as in the article
rule assemble_ilm:
    input:
        reads_r1 = lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "r1"),
        reads_r2 = lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "r2")
    wildcard_constraints:
        assembler = "|".join(ILM_ASSEMBLERS)
    params:
        outdir = "assembly/{assembler}",
        contig_name_format = "{sample}_{assembler}_%05d",
        assembler = "{assembler}"
    output:
        assembly = "assembly/{assembler}/{sample}_{assembler}.fna",
        log = "logs/assembly/{assembler}_{sample}.log" # ensure this is kept with shadow rule
    shadow: "shallow"
    threads: workflow.cores / 2 if (workflow.cores / 2 ) < 16 else 16 # shovill can't handle more
    conda: pathlib.Path(workflow.current_basedir) / "envs"/ "assembly_env.yaml"
    shell:
        """shovill --cpus {threads} \
            --R1 {input.reads_r1} \
            --R2 {input.reads_r2}\
            --force \
             --namefmt "{params.contig_name_format}" \
             --keepfiles \
             --assembler {params.assembler} \
             --trim \
            --outdir {params.outdir} \
            2> {output.log} || \
        touch "{params.outdir}/contigs.fa"
        mv "{params.outdir}/contigs.fa" "{output.assembly}" """

# ...but we may want to do some cleaning for the short read polishing anyway? TODO: cannibalize from RSYD-BASIC if so

# long read first: just stay with flye
rule trim_adapter:
    input:
        reads_ont = lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
                    input_data,
            "extra")
    output:
        reads_chopped = "reads/{sample}_trim.fastq.gz"
    conda: pathlib.Path(workflow.current_basedir) / "envs"/ "porechop_env.yaml"  # their version of porechop doesn't play nice with the remaining nanopore env
    shell:
        """porechop -i "{input.reads_ont}" -o "{output.reads_chopped}" """

#... but we'll do filtering with filtlong
rule filter_ont:
    input:
        reads_chopped = "reads/{sample}_trim.fastq.gz"
    output:
        reads_filtered = "reads/{sample}_filtered.fastq.gz"
    params:
        keep_percent = 95
    conda: pathlib.Path(workflow.current_basedir) / "envs"/ "nanopore_env.yaml"
    shell:
        """filtlong --keep_percent {params.keep_percent} "{input.reads_chopped}" | gzip > "{output.reads_filtered}" """

rule assemble_nanopore_only:
    input:
        filtered_fastq ="reads/{sample}_filtered.fastq.gz"
    output:
        nanopore_assembly = temp("assembly/flye_{sample}/assembly.fasta")
    params:
        out_dir = lambda wildcards, output: str(pathlib.Path(output.nanopore_assembly).parent)
    threads: workflow.cores / 2 if (workflow.cores / 2 ) < 16 else 16
    log: "logs/assembly/flye_{sample}.log"
    conda: pathlib.Path(workflow.current_basedir) / "envs"/ "nanopore_env.yaml"
    shadow: "shallow"
    shell:
        """
        flye --nano-raw "{input.filtered_fastq}" --out-dir "{params.out_dir}" \
            --threads {threads} &> "{log}" || touch "{output.nanopore_assembly}"
        """

rule medaka_consensus:
    input:
        filtered_fastq ="reads/{sample}_filtered.fastq.gz",
        flye_asm = "assembly/flye_{sample}/assembly.fasta"
    output:
        cleaned_asm = "assembly/flye_medaka/{sample}_flye_medaka.fna"
    threads: 16  # TODO: does medaka consensus actually claim four extra threads
    params:
        output_basedir = lambda wildcards, output: str(pathlib.Path(output.cleaned_asm).parent),
        snakemake_path = str(pathlib.Path(workflow.basedir).resolve()),
        medaka_model = "r941_min_hac_g507" # best we have
    log: "logs/assembly/medaka_{sample}.log"
    shadow: "shallow" # has to run as shadow so it doesn't trample on itself
    conda: pathlib.Path(workflow.current_basedir) / "envs" / "medaka_env.yaml"
    # only run this if we in fact have a sequence
    shell:
        """
        if [ ! -s {input.flye_asm} ]; then
            touch "{output.cleaned_asm}"
        else
            medaka_consensus -i "{input.filtered_fastq}" \
            -d "{input.flye_asm}" \
            -o "{params.output_basedir}" \
            -m {params.medaka_model} \
            -t {threads} -f > {log}
            mv "{params.output_basedir}/consensus.fasta" "{output.cleaned_asm}"
        fi
        """

rule index_asm:
    input:
        cleaned_asm = "assembly/flye_medaka/{sample}_flye_medaka.fna"
    output:
        asm_index="assembly/flye_medaka/{sample}_flye_medaka.fai"
    params:
        asm_prefix = lambda wildcards, output: str(pathlib.Path(output.asm_index).with_suffix(""))
    conda: pathlib.Path(workflow.current_basedir) / "envs" / "polish_env.yaml"
    shell:
        """bwa index "{input.cleaned_asm}" -p {params.asm_prefix} """


rule map_ilm_to_flye_asm:
    input:
        ilm_reads = lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
        input_data,
        wildcards.read_dir),
        cleaned_asm ="assembly/flye_medaka/{sample}_flye_medaka.fna",
        asm_index="assembly/flye_medaka/{sample}_flye_medaka.fai"
    output:
            alignment = "assembly/flye/{sample}_alignments_{read_dir}.sam"
    conda: pathlib.Path(workflow.current_basedir) / "envs"/"polish_env.yaml"
    threads: 4
    shell:
        """bwa mem -t {threads} -a "{input.cleaned_asm}" "{input.ilm_reads}" > "{output.alignment}" """

rule asm_polypolish:
    input:
        align_r1 = "assembly/flye/{sample}_alignments_r1.sam",
        align_r2 = "assembly/flye/{sample}_alignments_r2.sam",
        cleaned_asm = "assembly/flye_medaka/{sample}_flye_medaka.fna"
    output:
        polished = "assembly/lf_polypolish/{sample}_lf_polypolish.fna"
    conda: pathlib.Path(workflow.current_basedir) / "envs"/"polish_env.yaml"
    shell:
        """polypolish "{input.cleaned_asm}" "{input.align_r1}" "{input.align_r2}" > "{output.polished}" """



rule assemble_short_read_first:
    input:
        reads_r1=lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "r1"),
        reads_r2=lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "r2"),
        reads_ont = "reads/{sample}_filtered.fastq.gz"
    params:
        outdir = "assembly/sf_unicycler"
    output:
        outfile = "assembly/sf_unicycler/{sample}_sf_unicycler.fna"
    threads: workflow.cores / 2
    conda: pathlib.Path(workflow.current_basedir) /"envs"/ "assembly_env.yaml"
    log: "logs/assembly/sf_unicycler_{sample}.log"
    shadow: "shallow"  # need to have this because unicycler doesn't do output prefixes...
    shell:
        """
        unicycler -1 "{input.reads_r1}" -2 "{input.reads_r2}" --long "{input.reads_ont}" --threads {threads} \
         -o {params.outdir}
        mv "{params.outdir}/assembly.fasta" "{output.outfile}"
        mv "{params.outdir}/unicycler.log" "{log}"
        """



# they presumably have done this on raw reads, there's not much we can do there - TODO no they don't, we might try that too?
rule run_seqsero2:
    input:
        assembly  = "assembly/{assembler}/{sample}_{assembler}.fna"
    output:
        seqsero_report = "typing/seqsero2/{sample}_{assembler}_result.txt"
    shadow: "full"
    wildcard_constraints:
        assembler = "|".join(ALL_ASSEMBLERS)
    params:
        # input type: 2 is paired-end reads, 4 is assembly - do we want to tweak this?
        input_type = "4",
        workflow = "k",
        output_dir = "typing/seqsero2",
        placeholder_data = r"Sample name\tPredicted identification\tPredicted serotype\n" \
                           "{sample}" \
                           r"\t-\t-"
    threads: workflow.cores / 4 if (workflow.cores / 4) < 16 else 16
    log: "logs/serotyping/{assembler}_{sample}_seqsero2.log"
    conda: pathlib.Path(workflow.current_basedir) /"envs"/ "seqsero2_env.yaml"
    resources:
        mem_mb = 500  # todo: better estimate
    shell:
        """
        SeqSero2_package.py -t {params.input_type} -m {params.workflow} -p {threads} \
        -i {input.assembly} -d {params.output_dir} &> "{log}" || printf "{params.placeholder_data}" > "{params.output_dir}/SeqSero_result.txt"
        mv "{params.output_dir}/SeqSero_result.txt" "{output.seqsero_report}"
        """


rule run_sistr:
    input:
        assembly="assembly/{assembler}/{sample}_{assembler}.fna"
    output:
        sistr_report ="typing/sistr/{sample}_{assembler}_result.txt"
    threads: workflow.cores / 4 if (workflow.cores / 4) < 16 else 16
    wildcard_constraints:
        assembler="|".join(ALL_ASSEMBLERS)
    params:
        output_format = "tab"
    conda: pathlib.Path(workflow.current_basedir) / "envs"/"sistr_env.yaml"
    log: "logs/serotyping/{assembler}_{sample}_sistr.log"
    shell:
        """
        sistr "{input.assembly}" -f {params.output_format} -o "{output.sistr_report}" -t {threads}  &> "{log}"
        """