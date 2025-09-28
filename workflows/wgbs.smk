
# 1. Load the central configuration file
configfile: "config/wgbs_config.yaml"

# 2. Get a list of all sample names from the config file.
SAMPLES = list(config["samples"].keys())

# 3. Define the final desired output of the entire workflow.
rule all:
    input:
         expand("results/wgbs/qc/{sample}.done", sample=config["samples"])
        
        # The trimmed FASTQ files
        #expand("results/trimmed_reads/{sample}_R1_val_1.fq.gz", sample=SAMPLES),
        
        # QC reports for the new trimmed reads
        #expand("results/trimmed_qc/{sample}_R1_val_1_fastqc.html", sample=SAMPLES)

# 4. Include the modular rule files
include: "rules/qc.smk"
#include: "rules/trimming_and_qc.smk"
