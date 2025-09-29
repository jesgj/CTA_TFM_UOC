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

# Import all rules from wgbs module
use rule * from wgbs as wgbs_*

# Conditional rule_all must be defined at top level
if config["pipeline"] == "rnaseq":
    use rule all from rnaseq
elif config["pipeline"] == "wgbs":
    use rule all from wgbs as wgbs_all

    # Make it the default
    rule all:
        input: rules.wgbs_all.input
else:
    rule all:
        input: "ERROR: Write 'rnaseq' or 'wgbs' in config.yaml"
