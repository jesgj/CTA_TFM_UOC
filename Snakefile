# Load config file
configfile: "config/config.yaml"

# Import modules for different pipelines
module rnaseq:
    snakefile: "workflows/rnaseq.smk"
    config: config["rnaseq"]
    #prefix: "rnaseq_"

module wgbs:
    snakefile: "workflows/wgbs.smk"
    config: config["wgbs"]
    #prefix: "wgbs_"

# Conditional rule_all must be defined at top level, not inside a rule
if config["pipeline"] == "rnaseq":
    # Delegate 'all' rule to the module
    use rule all from rnaseq
elif config["pipeline"] == "wgbs":
    use rule all from wgbs
else:
    rule all:
        input: "ERROR: Write 'rnaseq' or 'wgbs' in config.yaml"
