# WGBS module workflow
import os
from snakemake.io import glob_wildcards

RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config.get("qc_dir", "results/wgbs/qc")

# Ensure QC directory exists
os.makedirs(QC_DIR, exist_ok=True)

# Detect sample IDs - capture everything before _R1
SAMPLES, = glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1_001.fastq.gz"))

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "{sample}_R1_001_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(QC_DIR, "{sample}_R2_001_fastqc.html"), sample=SAMPLES)