rule fastqc:
    input:
        R1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
        R2=lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        html_R1="results/{pipeline}/qc/{sample}_R1_fastqc.html",
        html_R2="results/{pipeline}/qc/{sample}_R2_fastqc.html",
        zip_R1="results/{pipeline}/qc/{sample}_R1_fastqc.zip",
        zip_R2="results/{pipeline}/qc/{sample}_R2_fastqc.zip",
        flag="results/{pipeline}/qc/{sample}.done"
    params:
         outdir = "results/{pipeline}/qc"
    log:
        "logs/fastqc_raw/{pipeline}{sample}.log"
    threads: 4
    shell:
        "pixi run fastqc {input.R1} {input.R2} --outdir {params.outdir} > {log} 2>&1"
        "&& touch {output.flag}"

