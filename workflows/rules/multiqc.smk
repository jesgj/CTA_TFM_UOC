# workflows/rules/multiqc.smk
import os

# This rule assumes the following config variables are set by the parent workflow:
# - pipeline_name: e.g., "chip_cr", "rnaseq", "wgbs"
# - multiqc_input_files: A list of all files that should be generated before multiqc runs.
# - multiqc_results_dir: The top-level results directory for the pipeline.

_MULTIQC_RESULTS_DIR = config["multiqc_results_dir"]
_PIPELINE_NAME = config["pipeline_name"]

# Construct the paths to scan. MultiQC is efficient at finding relevant files.
_ANALYSIS_DIRS = f"results/{_PIPELINE_NAME} logs/{_PIPELINE_NAME}"


rule multiqc:
    input:
        config["multiqc_input_files"]
    output:
        report = os.path.join(_MULTIQC_RESULTS_DIR, "multiqc_report.html")
    params:
        analysis_dirs = _ANALYSIS_DIRS,
        output_dir = _MULTIQC_RESULTS_DIR # Pass it to params
    log:
        os.path.join("logs", _PIPELINE_NAME, "multiqc.log")
    shell:
        "pixi run multiqc {params.analysis_dirs} -o {params.output_dir} --filename multiqc_report.html > {log}.out 2> {log}.err"
