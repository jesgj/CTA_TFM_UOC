
# Load config file
configfile: "config/config.yaml"

# Import modules for different pipelines
module rnaseq:
    snakefile: "workflows/rnaseq.smk"
    config: config["rnaseq"]
    prefix: "rnaseq_"

module wgbs:
    snakefile: "workflows/wgbs.smk"
    config: config["wgbs"]
    prefix: "wgbs_"

rule all:
    if config["pipeline"] == "rnaseq":
        use rule all from rnaseq
    elif config["pipeline"] == "wgbs":
        use rule all from wgbs
    else:
        input: "ERROR: Write 'rnaseq' o 'wgbs' in config.yaml"
