# WGBS module workflow
import os
import glob

RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config.get("qc_dir", "results/wgbs/qc")

# Make these available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR

# Ensure QC directory exists
os.makedirs(QC_DIR, exist_ok=True)

# Glob all fastq files to determine the final output for the 'all' rule.
ALL_FASTQ_FILES = glob.glob(os.path.join(RAW_DIR, '*.fastq.gz'))
BASENAMES = [os.path.basename(f).replace('.fastq.gz', '') for f in ALL_FASTQ_FILES]

print(f"Found {len(ALL_FASTQ_FILES)} fastq files.")

include: "rules/qc.smk"

rule all:
    input:
        expand(os.path.join(QC_DIR, "{basename}_fastqc.html"), basename=BASENAMES)
