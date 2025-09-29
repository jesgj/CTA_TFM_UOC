# WGBS module workflow
import os
from snakemake.io import glob_wildcards

RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config.get("qc_dir", "results/qc")  # fallback if not defined

# Ensure QC directory exists
os.makedirs(QC_DIR, exist_ok=True)

# Detect sample IDs from FASTQs (handles _R1/_R2 with optional _001)
SAMPLES, = glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1{extra}.fastq.gz"))

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(QC_DIR, "{sample}_R2_fastqc.html"), sample=SAMPLES)