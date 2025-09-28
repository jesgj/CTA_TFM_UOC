# WGBS module workflow
configfile: "config.yaml"

RAW_DIR = config["raw_fastqs_dir"]
SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/qc")  # fallback if not defined

include: "rules/qc.smk"

rule all: 
    input: 
        # Run QC on all samples 
        expand( f"{QC_DIR}/{{sample}}_R1_fastqc.html", sample=SAMPLES.keys() ), expand( f"{QC_DIR}/{{sample}}_R2_fastqc.html", sample=SAMPLES.keys() )