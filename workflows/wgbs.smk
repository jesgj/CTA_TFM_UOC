# WGBS module workflow
import os
import glob
from snakemake.io import glob_wildcards

# --- CONFIGURATION ---
# Define directories from config
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR

# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)


# --- MODULE INCLUSION ---

# 1. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"


# --- FINAL TARGETS ---

# Detect samples for final report generation
SAMPLES, = glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1_001.fastq.gz"))

# Glob all raw fastq files to determine the final output for the 'qc' part of the 'all' rule.
ALL_RAW_FASTQ_FILES = glob.glob(os.path.join(RAW_DIR, '*.fastq.gz'))
RAW_BASENAMES = [os.path.basename(f).replace('.fastq.gz', '') for f in ALL_RAW_FASTQ_FILES]


# Final target rule for the wgbs workflow
rule all:
    input:
        # 1. FastQC reports for raw files
        expand(os.path.join(QC_DIR, "{basename}_fastqc.html"), basename=RAW_BASENAMES),

        # 2. FastQC reports for trimmed files
        expand(os.path.join(QC_TRIMMED_DIR, "{sample}_{read}.trimmed_fastqc.html"),
               sample=SAMPLES,
               read=["R1_001", "R2_001"])