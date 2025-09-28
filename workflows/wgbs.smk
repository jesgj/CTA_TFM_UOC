# WGBS module workflow
RAW_DIR = config["raw_fastqs_dir"]
SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/qc")  # fallback if not defined

rule fastqc:
    input:
        R1=lambda wc: f"{RAW_DIR}/{SAMPLES[wc.sample]['R1']}",
        R2=lambda wc: f"{RAW_DIR}/{SAMPLES[wc.sample]['R2']}"
    output:
        html_R1=f"{QC_DIR}/{{sample}}_R1_fastqc.html",
        html_R2=f"{QC_DIR}/{{sample}}_R2_fastqc.html",
        zip_R1=f"{QC_DIR}/{{sample}}_R1_fastqc.zip",
        zip_R2=f"{QC_DIR}/{{sample}}_R2_fastqc.zip"
    threads: 2
    shell:
        "fastqc -o {QC_DIR} -t {threads} {input.R1} {input.R2}"