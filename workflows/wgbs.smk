# WGBS module workflow
import os
from snakemake.io import glob_wildcards

RAW_DIR = config["raw_fastqs_dir"]
#SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/qc")  # fallback if not defined

# Ensure QC directory exists
os.makedirs(QC_DIR, exist_ok=True)

# Detect all R1 FastQC HTML files dynamically
R1_files, = glob_wildcards(os.path.join(QC_DIR, "/*R1*fastqc.html"))
R2_files, = glob_wildcards(os.path.join(QC_DIR, "/*R2*fastqc.html"))

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "/{file}_fastqc.html"), file=R1_files),
        expand(os.path.join(QC_DIR, "/{file}_fastqc.html"), file=R2_files),