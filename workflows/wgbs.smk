# WGBS module workflow
import os
from snakemake.io import glob_wildcards

RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config.get("qc_dir", "results/qc")

# Ensure QC directory exists
os.makedirs(QC_DIR, exist_ok=True)

# Detect sample IDs from FASTQs (handles _R1/_R2 with optional _001)
# This captures the base sample name, stripping R1/R2 and suffix
SAMPLES_RAW, EXTRAS = glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1{extra}.fastq.gz"))

# Create a dictionary mapping sample to its extra suffix
SAMPLE_EXTRAS = dict(zip(SAMPLES_RAW, EXTRAS))

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "{sample}_R1{extra}_fastqc.html"), 
               zip, sample=SAMPLES_RAW, extra=EXTRAS),
        expand(os.path.join(QC_DIR, "{sample}_R2{extra}_fastqc.html"), 
               zip, sample=SAMPLES_RAW, extra=EXTRAS)