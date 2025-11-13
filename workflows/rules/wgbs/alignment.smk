# workflows/rules/alignment.smk
import os
from snakemake.io import touch

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
REF_DIR = os.path.dirname(REF_GENOME)
TRIMMED_DIR = config["trimmed_dir"]
ALIGN_DIR = config["alignment_dir"]
DEDUP_DIR = config["dedup_dir"]
BISMARK_EXTRA_ARGS = config.get("bismark", {}).get("extra_args", "")


# --- RULES ---

rule bismark_genome_preparation:
    """
    Creates a Bismark index for the reference genome if it doesn't exist.
    """
    input:
        REF_GENOME
    output:
        touch(os.path.join(REF_DIR, "Bisulfite_Genome", "bismark_index_created.OK"))
    params:
        ref_dir = REF_DIR
    threads: 1
    log:
        "logs/bismark_genome_preparation.log"
    shell:
        """
        pixi run bismark_genome_preparation --bowtie2 --verbose {params.ref_dir} > {log}.out 2> {log}.err
        """

rule bismark_alignment:
    """
    Aligns trimmed paired-end reads to the bisulfite genome using Bismark.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index = os.path.join(REF_DIR, "Bisulfite_Genome", "bismark_index_created.OK")
    output:
        bam = os.path.join(ALIGN_DIR, "{sample}_pe.bam"),
        report = os.path.join(ALIGN_DIR, "{sample}_pe_report.txt")
    params:
        ref_dir = REF_DIR,
        parallel = 8,
        align_dir = ALIGN_DIR,
        extra = BISMARK_EXTRA_ARGS
    threads: 64
    log:
        os.path.join("logs", "bismark", "{sample}.log")
    shell:
        """
        pixi run bismark --bowtie2 -p {params.parallel} {params.extra} \
        --genome {params.ref_dir} \
        -1 {input.r1} -2 {input.r2} \
        -o {params.align_dir} --basename {wildcards.sample} > {log}.out 2> {log}.err
        """

rule deduplicate_bismark:
    """
    Deduplicates a Bismark BAM file and moves it to the deduplication directory.
    """
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}_pe.bam")
    output:
        dedup_bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam"),
        report = os.path.join(DEDUP_DIR, "{sample}_pe.deduplication_report.txt")
    params:
        # The tool automatically creates this output file in the input directory
        auto_output_bam = lambda wildcards: os.path.join(ALIGN_DIR, f"{wildcards.sample}_pe.deduplicated.bam"),
        auto_output_report = lambda wildcards: os.path.join(ALIGN_DIR, f"{wildcards.sample}_pe.deduplication_report.txt")
    log:
        os.path.join("logs", "deduplicate_bismark", "{sample}.log")
    shell:
        """
        pixi run deduplicate_bismark -p --bam {input.bam} 2> {params.auto_output_report} 1> {log}
        mv {params.auto_output_bam} {output.dedup_bam}
        mv {params.auto_output_report} {output.report}
        """
