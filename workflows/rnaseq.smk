# workflows/rnaseq.smk
import os
import re
from collections import defaultdict

# --- CONFIGURATION ---
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
# Kallisto
KALLISTO_OUTPUT_DIR = config["kallisto_output_dir"]
TRANSCRIPTOME_FASTA = config["transcriptome_fasta"]
KALLISTO_INDEX = config["kallisto_index"]
# HISAT2
REF_GENOME = config["ref_genome"]
HISAT2_INDEX_DIR = config["hisat2_index_dir"]
ALIGNMENT_DIR = config["alignment_dir"]


# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["kallisto_output_dir"] = KALLISTO_OUTPUT_DIR
config["transcriptome_fasta"] = TRANSCRIPTOME_FASTA
config["kallisto_index"] = KALLISTO_INDEX
config["ref_genome"] = REF_GENOME
config["hisat2_index_dir"] = HISAT2_INDEX_DIR
config["alignment_dir"] = ALIGNMENT_DIR

# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(KALLISTO_OUTPUT_DIR, exist_ok=True)
os.makedirs(ALIGNMENT_DIR, exist_ok=True)
os.makedirs(os.path.join("logs", "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join("logs", "fastp"), exist_ok=True)
os.makedirs(os.path.join("logs", "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join("logs", "kallisto_quant"), exist_ok=True)
os.makedirs(os.path.join("logs", "hisat2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", "hisat2_build"), exist_ok=True)


# --- SAMPLE DISCOVERY ---
SAMPLES_INFO = config.get("samples_info", {})
SAMPLES = list(SAMPLES_INFO.keys())
config['samples_info'] = SAMPLES_INFO


# --- MODULE INCLUSION ---

# 1. QC on raw files (using the generic rule)
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files (using the generic rule)
include: "rules/trimming_and_qc.smk"

# 3. Pseudo-alignment
include: "rules/rnaseq/pseudoalignment.smk"

# 4. Alignment
include: "rules/rnaseq/alignment.smk"


# --- FINAL TARGETS ---

rule all:
    input:
        # Raw QC reports
        expand(os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"), sample=SAMPLES),
        # Trimmed QC reports
        expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.html"), sample=SAMPLES),
        # Kallisto quantification files
        expand(os.path.join(KALLISTO_OUTPUT_DIR, "{sample}", "abundance.tsv"), sample=SAMPLES),
        # HISAT2 alignment files
        expand(os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"), sample=SAMPLES)