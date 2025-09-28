RAW_DIR = config.get("raw_fastqs_dir")
SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/wgbs/qc")  # fallback default

rule fastqc:
    input:
        R1 = lambda wildcards: f"{RAW_DIR}/{SAMPLES[wildcards.sample]['R1']}",
        R2 = lambda wildcards: f"{RAW_DIR}/{SAMPLES[wildcards.sample]['R2']}"
    output:
        #html_R1 = QC_DIR + "/{sample}_R1_fastqc.html",
        #html_R2 = QC_DIR + "/{sample}_R2_fastqc.html",
        #zip_R1  = QC_DIR + "/{sample}_R1_fastqc.zip",
        #zip_R2  = QC_DIR + "/{sample}_R2_fastqc.zip"
        html_R1="results/wgbs/qc/{sample}_R1_fastqc.html",
        html_R2="results/wgbs/qc/{sample}_R2_fastqc.html",
        zip_R1="results/wgbs/qc/{sample}_R1_fastqc.zip",
        zip_R2="results/wgbs/qc/{sample}_R2_fastqc.zip"
    threads: 2
    shell:
        "pixi run fastqc -o {QC_DIR} -t {threads} {input.R1} {input.R2}"

