# WGBS module workflow
import os
import re

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

# --- SAMPLE DISCOVERY ---
def discover_samples_and_reads(raw_dir):
    """
    Scans the raw directory for paired-end fastq files and returns a dictionary
    mapping sample names to their R1 and R2 file paths.
    Handles both _R1.fq.gz and _R1_001.fq.gz naming conventions.
    """
    samples_info = {}
    if not os.path.exists(raw_dir):
        return samples_info

    for f in os.listdir(raw_dir):
        # The regex now looks for .fq.gz and captures the sample name and read identifier part
        match = re.match(r'(.+?)(_R1|_R1_001)\.fq\.gz', f)
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

# Make sample info available to included rules
config['samples_info'] = SAMPLES_INFO


# --- MODULE INCLUSION ---

# 1. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"


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
               read=["R1", "R2"])
