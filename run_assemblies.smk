import pathlib

import pandas as pd

import snake_helpers

# read file overview: sample name, path to r1, path to r2, extra (ONT reads)
input_data = pd.read_csv(config["infile"], sep="\t")

ILM_ASSEMBLERS = {"skesa"}
HYBRID_ASSEMBLERS = {"sf_unicycler"}  # , "lf_polypolish"}
ALL_ASSEMBLERS = ILM_ASSEMBLERS.union(HYBRID_ASSEMBLERS)
SEROTYPERS = {"seqsero2"} #  {"seqsero", "seqsero2", "sistr"} TODO add again

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

rule assemble_short_read_first:
    input:
        reads_r1=lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "r1"),
        reads_r2=lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "r2"),
        reads_ont = lambda wildcards: snake_helpers.get_reads_from_overview(wildcards.sample,
            input_data,
            "extra")
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

# they presumably have done this on raw reads, there's not much we can do there
rule run_seqsero2:
    input:
        assembly  = "assembly/{assembler}/{sample}_{assembler}.fna"
    output:
        seqsero_report = "typing/SeqSero2/{sample}_{assembler}_result.txt"
    shadow: "full"
    conda: "typing_env"
    params:
        # input type: 2 is paired-end reads, 4 is assembly - do we want to tweak this?
        input_type = "4",
        workflow = "k",
        output_dir = "typing/SeqSero2",
        placeholder_data = r"Sample name\tPredicted identification\tPredicted serotype\n" \
                           "{sample}" \
                           r"\t-\t-"
    threads: workflow.cores / 4 if (workflow.cores / 4) < 16 else 16
    conda: pathlib.Path(workflow.current_basedir) / "seqsero2_env.yaml"
    resources:
        mem_mb = 500  # todo: better estimate
    shell:
        """
        SeqSero2_package.py -t {params.input_type} -m {params.workflow} -p {threads} \
        -i {input.assembly} -d {params.output_dir} || printf "{params.placeholder_data}" > "{params.output_dir}/SeqSero_result.txt"
        mv "{params.output_dir}/SeqSero_result.txt" "{output.seqsero_report}"
        """
