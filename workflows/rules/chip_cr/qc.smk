# workflows/rules/chip_cr/qc.smk
import os

SAMPLES_INFO = config['samples_info']
QC_DIR = config["qc_dir"]

rule fastqc_raw_pe:
    """
    Runs FastQC on raw paired-end reads for a single sample.
    This rule is selected for PE samples because rule 'all' requests
    the specific PE output files.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1'],
        r2 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R2']
    output:
        html_r1 = os.path.join(QC_DIR, "{sample}_R1_fastqc.html"),
        html_r2 = os.path.join(QC_DIR, "{sample}_R2_fastqc.html"),
        zip_r1 = os.path.join(QC_DIR, "{sample}_R1_fastqc.zip"),
        zip_r2 = os.path.join(QC_DIR, "{sample}_R2_fastqc.zip")
    params:
        outdir = QC_DIR
    threads: 2
    log:
        os.path.join("logs", "chipseq_cutrun", "fastqc_raw", "{sample}_pe.log")
    shell:
        "pixi run fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} &> {log}"

rule fastqc_raw_se:
    """
    Runs FastQC on raw single-end reads for a single sample.
    This rule is selected for SE samples because rule 'all' requests
    the specific SE output files.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1']
    output:
        html = os.path.join(QC_DIR, "{sample}_fastqc.html"),
        zip = os.path.join(QC_DIR, "{sample}_fastqc.zip")
    params:
        outdir = QC_DIR
    threads: 1
    log:
        os.path.join("logs", "chipseq_cutrun", "fastqc_raw", "{sample}_se.log")
    shell:
        "pixi run fastqc -o {params.outdir} -t {threads} {input.r1} &> {log}"