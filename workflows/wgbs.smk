# FILE: workflows/wgbs.smk
# Variables from config
RAW_DIR = config.get("raw_fastqs_dir")
SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/wgbs/qc")  # Added fallback default

print(f"RAW_DIR: {RAW_DIR}")
print(f"SAMPLES: {SAMPLES}")
print(f"QC_DIR: {QC_DIR}")

# Include shared QC rules
include: "rules/qc.smk"

# This rule defines the final output of the pipeline
rule all:
    input:
        expand("results/wgbs/qc/{sample}_R1_fastqc.html", sample=SAMPLES.keys()),
        expand("results/wgbs/qc/{sample}_R2_fastqc.html", sample=SAMPLES.keys())
    shell:
         "echo DONUTS"
