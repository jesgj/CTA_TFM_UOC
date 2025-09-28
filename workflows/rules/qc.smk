RAW_DIR = config.get("raw_fastqs_dir")
SAMPLES = config["samples"]
QC_DIR = config.get("qc_dir", "results/wgbs/qc")  # fallback default

rule fastqc:
    input:
        R1 = lambda wc: f"{RAW_DIR}/{SAMPLES[wc.sample]['R1']}",
        R2 = lambda wc: f"{RAW_DIR}/{SAMPLES[wc.sample]['R2']}"
    output:
        html_R1 = lambda wc: f"{QC_DIR}/{Path(SAMPLES[wc.sample]['R1']).stem}_fastqc.html",
        html_R2 = lambda wc: f"{QC_DIR}/{Path(SAMPLES[wc.sample]['R2']).stem}_fastqc.html",
        zip_R1  = lambda wc: f"{QC_DIR}/{Path(SAMPLES[wc.sample]['R1']).stem}_fastqc.zip",
        zip_R2  = lambda wc: f"{QC_DIR}/{Path(SAMPLES[wc.sample]['R2']).stem}_fastqc.zip"
    threads: 4
    shell:
        "pixi run fastqc -o {QC_DIR} -t {threads} {input.R1} {input.R2}"

