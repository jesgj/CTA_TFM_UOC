import os

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