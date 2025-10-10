# rules/qc.smk - Generic QC rules (refactored for per-sample execution)
import os

# This module now expects SAMPLES_INFO to be in the config
SAMPLES_INFO = config['samples_info']
QC_DIR = config["qc_dir"]

rule fastqc_raw:
    """
    Runs FastQC on raw paired-end reads for a single sample.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1'],
        r2 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R2']
    output:
        html_r1 = os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"),
        html_r2 = os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"),
        zip_r1 = os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.zip"),
        zip_r2 = os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.zip")
    params:
        outdir = QC_DIR
    threads: 4 # fastqc can use multiple threads for multiple files
    log:
        os.path.join("logs", "fastqc_raw", "{sample}.log")
    shell:
        """
        pixi run fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} &> {log}
        """