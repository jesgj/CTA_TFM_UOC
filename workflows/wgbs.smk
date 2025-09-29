# WGBS module workflow
import os
from snakemake.io import glob_wildcards

# Get absolute paths to verify
WORKFLOW_DIR = os.getcwd()
print(f"Working directory: {WORKFLOW_DIR}")

RAW_DIR = config["wgbs"]["raw_fastqs_dir"]
QC_DIR = config["wgbs"]["qc_dir"]

print(f"RAW_DIR (relative): {RAW_DIR}")
print(f"RAW_DIR (absolute): {os.path.abspath(RAW_DIR)}")
print(f"RAW_DIR exists: {os.path.exists(RAW_DIR)}")
print(f"QC_DIR (relative): {QC_DIR}")
print(f"QC_DIR (absolute): {os.path.abspath(QC_DIR)}")

# List files in RAW_DIR to verify
if os.path.exists(RAW_DIR):
    files = os.listdir(RAW_DIR)
    print(f"Files in RAW_DIR: {files[:5]}...")  # Show first 5

# Ensure QC directory exists
os.makedirs(QC_DIR, exist_ok=True)

# Detect sample IDs
SAMPLES, = glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1_001.fastq.gz"))

print(f"Found {len(SAMPLES)} samples: {SAMPLES}")

# Test one sample path construction
if SAMPLES:
    test_sample = SAMPLES[0]
    test_input = os.path.join(RAW_DIR, f"{test_sample}_R1_001.fastq.gz")
    test_output = os.path.join(QC_DIR, f"{test_sample}_R1_001_fastqc.html")
    print(f"\nTest paths for sample '{test_sample}':")
    print(f"  Input:  {test_input}")
    print(f"  Input exists: {os.path.exists(test_input)}")
    print(f"  Output: {test_output}")

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "{sample}_R1_001_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(QC_DIR, "{sample}_R2_001_fastqc.html"), sample=SAMPLES)