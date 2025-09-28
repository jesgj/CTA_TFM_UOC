# Load the configuration file passed by the user
configfile: workflow.configfile

PIPELINE = config.get("pipeline", None)

if PIPELINE == "wgbs":
    include: "workflows/wgbs.smk"
elif PIPELINE == "rnaseq":
    include: "workflows/rnaseq.smk"
else:
    print(f"ERROR: Unknown pipeline '{PIPELINE}'")
    exit(1)