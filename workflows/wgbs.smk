# WGBS module workflow
import os
import re

# --- CONFIGURATION ---
# Define directories from config
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
ALIGN_DIR = config["alignment_dir"]
DEDUP_DIR = config["dedup_dir"]
REF_GENOME = config["ref_genome"]

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["alignment_dir"] = ALIGN_DIR
config["dedup_dir"] = DEDUP_DIR
config["ref_genome"] = REF_GENOME

# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(ALIGN_DIR, exist_ok=True)
os.makedirs(DEDUP_DIR, exist_ok=True)

# --- SAMPLE DISCOVERY ---
def discover_samples_and_reads(raw_dir):
    samples_info = {}
    if not os.path.exists(raw_dir):
        return samples_info
    for f in os.listdir(raw_dir):
        match = re.match(r'(.+?)(_R1|_R1_001)\.(fq|fastq)\.gz$', f)
        if match:
            sample_name = match.group(1)
            r1_part = match.group(2)
            r2_part = r1_part.replace('_R1', '_R2')
            r1_path = os.path.join(raw_dir, f)
            r2_path = os.path.join(raw_dir, f.replace(r1_part, r2_part))
            if os.path.exists(r2_path):
                samples_info[sample_name] = {'R1': r1_path, 'R2': r2_path}
    return samples_info

SAMPLES_INFO = discover_samples_and_reads(RAW_DIR)
SAMPLES = list(SAMPLES_INFO.keys())

config['samples_info'] = SAMPLES_INFO


# --- MODULE INCLUSION ---

# 1. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"

# 3. Alignment
include: "rules/alignment.smk"


# --- FINAL TARGETS ---

# Final target rule for the wgbs workflow
rule all:
    input:
        # 1. FastQC reports for raw files
        expand(os.path.join(QC_DIR, "{sample}_{read}_raw_fastqc.html"),
               sample=SAMPLES,
               read=["R1", "R2"]),

        # 2. FastQC reports for trimmed files
        expand(os.path.join(QC_TRIMMED_DIR, "{sample}_{read}_trimmed_fastqc.html"),
               sample=SAMPLES,
               read=["R1", "R2"]),
        
        # 3. Bismark deduplicated alignment files
        expand(os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam"), sample=SAMPLES)
