import os, re
from pathlib import Path

RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config.get("qc_dir", "results/wgbs/qc")  # fallback default

def build_map():
    fastqs = [f for f in os.listdir(RAW_DIR) if f.endswith((".fastq.gz", ".fastq"))]
    samples = {}
    for fq in fastqs:
        # match ..._R1.fastq.gz or ..._R1_001.fastq.gz (same for R2)
        m = re.match(r"(.+)_R([12])(?:_001)?\.f(ast)?q(\.gz)?$", fq)
        if not m:
            raise ValueError(f"Unrecognized FASTQ naming: {fq}")
        sample_id, read = m.group(1), m.group(2)
        if sample_id not in samples:
            samples[sample_id] = {}
        samples[sample_id][f"R{read}"] = fq
    return samples    

SAMPLES = build_map()

rule fastqc:
    input:
        R1 = lambda wc: os.path.join(RAW_DIR, SAMPLES[wc.sample]["R1"]),
        R2 = lambda wc: os.path.join(RAW_DIR, SAMPLES[wc.sample]["R2"])
    output:
        html_R1 = f"{QC_DIR}/{{sample}}_R1_fastqc.html",
        html_R2 = f"{QC_DIR}/{{sample}}_R2_fastqc.html",
        zip_R1  = f"{QC_DIR}/{{sample}}_R1_fastqc.zip",
        zip_R2  = f"{QC_DIR}/{{sample}}_R2_fastqc.zip"
    threads: 4
    shell:
        "pixi run fastqc -o {QC_DIR} -t {threads} {input.R1} {input.R2}"

