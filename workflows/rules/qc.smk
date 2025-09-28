
# RAW_DIR and SAMPLES come from the module's config
RAW_DIR = config.get("raw_fastqs_dir")
SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/qc")

rule fastqc:
    """
    Run FastQC on paired-end reads for all samples.
    """
    input:
        R1=lambda wildcards: f"{RAW_DIR}/{SAMPLES[wildcards.sample]['R1']}",
        R2=lambda wildcards: f"{RAW_DIR}/{SAMPLES[wildcards.sample]['R2']}"
    output:
        html_R1="results/qc/{sample}_R1_fastqc.html",
        html_R2="results/qc/{sample}_R2_fastqc.html",
        zip_R1="results/qc/{sample}_R1_fastqc.zip",
        zip_R2="results/qc/{sample}_R2_fastqc.zip"
    threads: 2
    shell:
        ""pixi run fastqc -o results/qc -t {threads} {input.R1} {input.R2}""