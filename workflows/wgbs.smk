# WGBS module workflow
import os
from snakemake.io import glob_wildcards

RAW_DIR = config["raw_fastqs_dir"]
#SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/qc")  # fallback if not defined

# Detect all R1 FastQC HTML files dynamically
R1_files, = glob_wildcards(os.path.join(QC_DIR, "*_R1_*_fastqc.html"))
R2_files, = glob_wildcards(os.path.join(QC_DIR, "*_R2_*_fastqc.html"))

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "{file}_fastqc.html"), file=R1_files),
        expand(os.path.join(QC_DIR, "{file}_fastqc.html"), file=R2_files),
        expand(os.path.join(QC_DIR, "{file}_fastqc.zip"),  file=R1_files),
        expand(os.path.join(QC_DIR, "{file}_fastqc.zip"),  file=R2_files)