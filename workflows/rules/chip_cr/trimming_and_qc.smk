# workflows/rules/chip_cr/trimming_and_qc.smk
import os

SAMPLES_INFO = config['samples_info']
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
FASTP_EXTRA_ARGS = config.get("fastp", {}).get("extra_args", "")

# --- Trimming ---

rule fastp_trim_pe:
    """
    Trims paired-end reads using fastp.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1'],
        r2 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R2']
    output:
        trimmed_r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        trimmed_r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        html = os.path.join(TRIMMED_DIR, "{sample}_pe.fastp.html"),
        json = os.path.join(TRIMMED_DIR, "{sample}_pe.fastp.json")
    params:
        extra = FASTP_EXTRA_ARGS
    threads: 4
    log:
        os.path.join("logs", "chipseq_cutrun", "fastp", "{sample}_pe.log")
    shell:
        """
        pixi run fastp -i {input.r1} -I {input.r2} \
        -o {output.trimmed_r1} -O {output.trimmed_r2} \
        -h {output.html} -j {output.json} \
        {params.extra} \
        -t {threads} &> {log}
        """

rule fastp_trim_se:
    """
    Trims single-end reads using fastp.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1']
    output:
        trimmed_r1 = os.path.join(TRIMMED_DIR, "{sample}_SE.trimmed.fq.gz"),
        html = os.path.join(TRIMMED_DIR, "{sample}_se.fastp.html"),
        json = os.path.join(TRIMMED_DIR, "{sample}_se.fastp.json")
    params:
        extra = FASTP_EXTRA_ARGS
    threads: 4
    log:
        os.path.join("logs", "chipseq_cutrun", "fastp", "{sample}_se.log")
    shell:
        """
        pixi run fastp -i {input.r1} \
        -o {output.trimmed_r1} \
        -h {output.html} -j {output.json} \
        {params.extra} \
        -t {threads} &> {log}
        """

# --- QC on Trimmed ---

rule fastqc_trimmed_pe:
    """
    Runs FastQC on trimmed paired-end reads.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz")
    output:
        html_r1 = os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.html"),
        html_r2 = os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.html"),
        zip_r1 = os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.zip"),
        zip_r2 = os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.zip")
    params:
        outdir = QC_TRIMMED_DIR
    threads: 2
    log:
        os.path.join("logs", "chipseq_cutrun", "fastqc_trimmed", "{sample}_pe.log")
    shell:
        "pixi run fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} &> {log}"


rule fastqc_trimmed_se:
    """
    Runs FastQC on trimmed single-end reads.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_SE.trimmed.fq.gz")
    output:
        html = os.path.join(QC_TRIMMED_DIR, "{sample}_SE_trimmed_fastqc.html"),
        zip = os.path.join(QC_TRIMMED_DIR, "{sample}_SE_trimmed_fastqc.zip")
    params:
        outdir = QC_TRIMMED_DIR
    threads: 1
    log:
        os.path.join("logs", "chipseq_cutrun", "fastqc_trimmed", "{sample}_se.log")
    shell:
        "pixi run fastqc -o {params.outdir} -t {threads} {input.r1} &> {log}"
