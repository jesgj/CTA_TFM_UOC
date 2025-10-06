# workflows/rules/trimming_and_qc.smk
import os
from snakemake.io import glob_wildcards

# This module expects the following config variables to be set by the parent workflow:
# - raw_fastqs_dir: Path to the directory with raw fastq files.
# - trimmed_dir: Path to the directory where trimmed fastq files will be stored.
# - qc_trimmed_dir: Path to the directory for FastQC reports of trimmed files.
#-------------------------------------------------------------------------------------------------
# I have created the new module at workflows/rules/trimming_and_qc.smk.

#To use it, you will need to include it in a parent Snakefile (like workflows/wgbs.smk) and define the required config variables (raw_fastqs_dir, trimmed_dir, and qc_trimmed_dir) before the 
#include statement. You'll also need to update your rule all to request the final files you want to generate, for example the HTML reports from fastqc_trimmed.
#-------------------------------------------------------------------------------------------------
RAW_DIR = config["raw_fastqs_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]

# Detect samples from raw fastq files
SAMPLES, = glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1_001.fastq.gz"))

rule fastp_trim:
    """
    Trims paired-end reads using fastp.
    """
    input:
        r1 = os.path.join(RAW_DIR, "{sample}_R1_001.fastq.gz"),
        r2 = os.path.join(RAW_DIR, "{sample}_R2_001.fastq.gz")
    output:
        trimmed_r1 = os.path.join(TRIMMED_DIR, "{sample}_R1_001.trimmed.fastq.gz"),
        trimmed_r2 = os.path.join(TRIMMED_DIR, "{sample}_R2_001.trimmed.fastq.gz"),
        html = os.path.join(TRIMMED_DIR, "{sample}.fastp.html"),
        json = os.path.join(TRIMMED_DIR, "{sample}.fastp.json")
    threads: 4
    log:
        os.path.join("logs", "fastp", "{sample}.log")
    shell:
        """
        pixi run fastp -i {input.r1} -I {input.r2} \
        -o {output.trimmed_r1} -O {output.trimmed_r2} \
        -h {output.html} -j {output.json} \
        -t {threads} &> {log}
        """

rule fastqc_trimmed:
    """
    Runs FastQC on trimmed paired-end reads.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1_001.trimmed.fastq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2_001.trimmed.fastq.gz")
    output:
        html_r1 = os.path.join(QC_TRIMMED_DIR, "{sample}_R1_001.trimmed_fastqc.html"),
        html_r2 = os.path.join(QC_TRIMMED_DIR, "{sample}_R2_001.trimmed_fastqc.html"),
        zip_r1 = os.path.join(QC_TRIMMED_DIR, "{sample}_R1_001.trimmed_fastqc.zip"),
        zip_r2 = os.path.join(QC_TRIMMED_DIR, "{sample}_R2_001.trimmed_fastqc.zip")
    params:
        outdir = QC_TRIMMED_DIR
    threads: 4
    log:
        os.path.join("logs", "fastqc_trimmed", "{sample}.log")
    shell:
        """
        pixi run fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} &> {log}
        """

