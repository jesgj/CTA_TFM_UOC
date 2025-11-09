# workflows/rules/qc.smk
import os

# This module now expects SAMPLES_INFO to be in the config
SAMPLES_INFO = config['samples_info']
QC_DIR = config["qc_dir"]

rule fastqc_raw_pe:
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
    threads: 2
    log:
        os.path.join("logs", config["pipeline"], "fastqc_raw", "{sample}_pe.log")
    shell:
        """
        pixi run fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} &> {log}
        mv {params.outdir}/{wildcards.sample}_R1.fastqc.html {output.html_r1}
        mv {params.outdir}/{wildcards.sample}_R2.fastqc.html {output.html_r2}
        mv {params.outdir}/{wildcards.sample}_R1.fastqc.zip {output.zip_r1}
        mv {params.outdir}/{wildcards.sample}_R2.fastqc.zip {output.zip_r2}
        """

rule fastqc_raw_se:
    """
    Runs FastQC on raw single-end reads for a single sample.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1']
    output:
        html = os.path.join(QC_DIR, "{sample}_raw_fastqc.html"),
        zip = os.path.join(QC_DIR, "{sample}_raw_fastqc.zip")
    wildcard_constraints:
        sample=r"^(?!.*_R[12]$).*"
    params:
        outdir = QC_DIR
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "fastqc_raw", "{sample}_se.log")
    shell:
        """
        pixi run fastqc -o {params.outdir} -t {threads} {input.r1} &> {log}
        mv {params.outdir}/{wildcards.sample}.fastqc.html {output.html}
        mv {params.outdir}/{wildcards.sample}.fastqc.zip {output.zip}
        """
