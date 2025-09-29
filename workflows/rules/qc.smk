# rules/qc.smk - Generic QC rules
import os

# These should be set by the parent workflow before including this file
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]

print(f"RAW_DIR (relative): {RAW_DIR}")
print(f"RAW_DIR (absolute): {os.path.abspath(RAW_DIR)}")
print(f"RAW_DIR exists: {os.path.exists(RAW_DIR)}")
print(f"QC_DIR (relative): {QC_DIR}")
print(f"QC_DIR (absolute): {os.path.abspath(QC_DIR)}")

# List files in RAW_DIR to verify
if os.path.exists(RAW_DIR):
    files = os.listdir(RAW_DIR)
    print(f"Files in RAW_DIR: {files[:5]}...")  # Show first 5

rule fastqc:
    input:
        R1 = os.path.join(RAW_DIR, "{sample}_R1_001.fastq.gz"),
        R2 = os.path.join(RAW_DIR, "{sample}_R2_001.fastq.gz")
    output:
        html_R1 = os.path.join(QC_DIR, "{sample}_R1_001_fastqc.html"),
        html_R2 = os.path.join(QC_DIR, "{sample}_R2_001_fastqc.html"),
        zip_R1  = os.path.join(QC_DIR, "{sample}_R1_001_fastqc.zip"),
        zip_R2  = os.path.join(QC_DIR, "{sample}_R2_001_fastqc.zip")
    params:
        outdir = QC_DIR
    threads: 4
    shell:
        """
        pixi run fastqc -o {params.outdir} -t {threads} {input.R1} {input.R2}
        """