# rules/qc.smk - Generic QC rules
import os
import glob

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

# Use glob to find all fastq.gz files.
# This list is used to define the output files for the 'fastqc' rule.
ALL_FASTQ_FILES = glob.glob(os.path.join(RAW_DIR, '*.fastq.gz'))
BASENAMES = [os.path.basename(f).replace('.fastq.gz', '') for f in ALL_FASTQ_FILES]

rule fastqc:
    input:
        ALL_FASTQ_FILES
    output:
        htmls = expand(os.path.join(QC_DIR, '{basename}_fastqc.html'), basename=BASENAMES),
        zips  = expand(os.path.join(QC_DIR, '{basename}_fastqc.zip'), basename=BASENAMES)
    params:
        outdir = QC_DIR
    threads: 4
    shell:
        """
        pixi run fastqc -o {params.outdir} -t {threads} {input}
        """
