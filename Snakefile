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

module chip_cr:
    snakefile: "workflows/chip_cr.smk"
    config: config["chip_cr"]

# Conditional rule import and rule_all definition
if config["pipeline"] == "rnaseq":
    use rule * from rnaseq
    use rule all from rnaseq
elif config["pipeline"] == "wgbs":
    use rule * from wgbs as wgbs_*
    use rule all from wgbs as wgbs_all

    # Make it the default
    rule all:
        input: rules.wgbs_all.input
elif config["pipeline"] == "chip_cr":
    use rule * from chip_cr as chip_cr_*
    use rule all from chip_cr as chip_cr_all

    rule all:
        input: rules.chip_cr_all.input
else:
    rule all:
        input: "ERROR: Write 'rnaseq', 'wgbs', or 'chip_cr' in config.yaml"
