import pathlib

import pandas as pd

import snake_helpers

# read file overview: sample name, path to r1, path to r2, extra (ONT reads)
input_data = pd.read_csv(config["infile"], sep="\t")

ILM_ASSEMBLERS = {"skesa"}
HYBRID_ASSEMBLERS = {"sf_unicycler"}  # , "lf_polypolish"}
ALL_ASSEMBLERS = ILM_ASSEMBLERS.union(HYBRID_ASSEMBLERS)
SEROTYPERS = {"seqsero2", "sistr"}

rule all:
    input:
        all_serotyped = expand("typing/{serotyper}/{sample}_{assembler}_result.txt", # we need to move the result
        serotyper = SEROTYPERS, sample=input_data["sample"])

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
    shadow: "full"
    threads: workflow.cores / 2 if (workflow.cores / 2 ) < 16 else 16 # shovill can't handle more
    conda: pathlib.Path(workflow.current_basedir) / "assembly_env.yaml"
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
    conda: pathlib.Path(workflow.current_basedir) / "nanopore_env.yaml"
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
    conda: pathlib.Path(workflow.current_basedir) / "nanopore_env.yaml"
    shell:
        """filtlong --keep_percent {params.keep_percent} "{input.reads_chopped}" | gzip > "{output.reads_filtered}" """


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
    conda: pathlib.Path(workflow.current_basedir) / "assembly_env.yaml"
    log: "logs/assembly/sf_unicycler_{sample}.log"
    shadow: "full"  # need to have this because unicycler doesn't do output prefixes...
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
        seqsero_report = "typing/SeqSero2/{sample}_{assembler}_result.txt"
    shadow: "full"
    params:
        # input type: 2 is paired-end reads, 4 is assembly - do we want to tweak this?
        input_type = "4",
        workflow = "k",
        output_dir = "typing/SeqSero2",
        placeholder_data = r"Sample name\tPredicted identification\tPredicted serotype\n" \
                           "{sample}" \
                           r"\t-\t-"
    threads: workflow.cores / 4 if (workflow.cores / 4) < 16 else 16
    log: "logs/serotyping/{assembler}_{sample}_seqsero2.log"
    conda: pathlib.Path(workflow.current_basedir) / "seqsero2_env.yaml"
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
    params:
        output_format = "tab"
    conda: pathlib.Path(workflow.current_basedir) / "sistr_env.yaml"
    log: "logs/serotyping/{assembler}_{sample}_sistr.log"
    shell:
        """
        sistr_cmd "{input.assembly}" -f {params.output_format} -o "{output.sistr_report}" -t {threads}  &> "{log}"
        """