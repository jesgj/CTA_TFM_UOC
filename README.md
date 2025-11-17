# Modular Bioinformatics Pipeline for Next-Generation Sequencing Analysis

This project provides a comprehensive and reproducible framework for analyzing Next-Generation Sequencing (NGS) data using [Snakemake](https://snakemake.readthedocs.io/en/stable/). It is designed to be highly modular, allowing users to run different analysis workflows by simply changing a setting in a configuration file.

The project currently supports three main pipelines:
- **ChIP-seq / CUT&RUN (`chip_cr`)**: For analyzing protein-DNA interactions.
- **Whole-Genome Bisulfite Sequencing (`wgbs`)**: For analyzing DNA methylation.
- **RNA-seq (`rnaseq`)**: For transcriptomics and gene expression analysis.

---

## Project Overview

This project is a highly modular and reproducible bioinformatics pipeline for analyzing Next-Generation Sequencing (NGS) data. It is built using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system and utilizes [Pixi](https://pixi.sh/) for dependency and environment management.

The primary goal of this project is to provide a flexible framework for running different types of NGS analysis by simply changing a configuration setting. It currently supports three distinct pipelines:

*   **ChIP-seq / CUT&RUN (`chip_cr`)**: For analyzing protein-DNA interactions.
*   **Whole-Genome Bisulfite Sequencing (`wgbs`)**: For analyzing DNA methylation.
*   **RNA-seq (`rnaseq`)**: For transcriptomics and gene expression analysis.

### Key Features

*   **Modularity**: The pipeline is architected around a main `Snakefile` that acts as a controller, dynamically loading the selected workflow (`chip_cr`, `wgbs`, or `rnaseq`) from the `workflows/` directory.
*   **Configuration-driven**: All aspects of the pipeline, from sample definitions to tool-specific parameters, are controlled through a central configuration file, `config/config.yaml`.
*   **Reproducibility**: The project uses `pixi.toml` to define and lock the versions of all bioinformatics tools and software dependencies. This ensures that the analysis can be reproduced consistently across different environments.
*   **Automated Sample Discovery**: The pipeline can automatically discover input FASTQ files from a directory (as defined in `config.yaml`) or use a manually specified sample sheet. This is handled by helper functions in `src/utils.py`.
*   **Unified QC**: A centralized set of rules for quality control (`fastqc`) and adapter trimming (`fastp`) are shared across all pipelines, ensuring consistent data preprocessing.
*   **Comprehensive Reporting**: Each pipeline culminates in a [MultiQC](https://multiqc.info/) report, which aggregates results from all analysis steps into a single, interactive HTML file.

## Building and Running the Pipeline

The pipeline is executed using the `snakemake` command, which should be run within the Pixi environment to ensure all dependencies are available.

### 1. Setup

The project dependencies are managed by Pixi. To install them, you would typically run:

```bash
pixi install
```

### 2. Configuration

1.  **Select the pipeline**: Open `config/config.yaml` and set the `pipeline` variable to the desired workflow: `"chip_cr"`, `"wgbs"`, or `"rnaseq"`.
2.  **Configure paths and samples**:
    *   Update the `raw_fastqs_dir` to point to your input data.
    *   Either provide a `samples_info` block with a list of your samples and their FASTQ files or leave it empty to let the pipeline auto-discover samples from the `raw_fastqs_dir`.
    *   Adjust any other pipeline-specific parameters as needed (e.g., reference genome path, tool arguments).

### 3. Execution

To run the pipeline, use the `snakemake` command. It is recommended to run it through Pixi to ensure the correct environment is used.

**Example commands:**

*   **Dry-run (to see what tasks will be executed):**

    ```bash
    pixi run snakemake -- -n
    ```

*   **Execute the pipeline locally (using 8 cores):**

    ```bash
    pixi run snakemake -- -c8
    ```

*   **Generate a DAG (Directed Acyclic Graph) of the workflow:**

    ```bash
    pixi run snakemake --dag | dot -Tpng > dag.png
    ```

## Development Conventions

*   **Workflow Structure**: Each main pipeline (e.g., `chip_cr`) has a master snakefile in `workflows/` (e.g., `workflows/chip_cr.smk`). This file is responsible for loading configuration, discovering samples, and including the necessary rule files.
*   **Modular Rules**: The actual pipeline logic is encapsulated in modular rule files located in `workflows/rules/`.
    *   Generic rules that can be shared across pipelines (e.g., `qc.smk`, `multiqc.smk`) are in the top-level `rules/` directory.
    *   Pipeline-specific rules are organized into subdirectories (e.g., `workflows/rules/chip_cr/`).
*   **Configuration**: All configurable parameters should be in `config/config.yaml`. The Snakefile rules should read their parameters from the `config` object.
*   **Scripts**: Custom helper scripts and utility functions are placed in the `src/` directory.
*   **Logging**: All rules are configured to write `stdout` and `stderr` to log files in the `logs/` directory, with a subdirectory for each pipeline and rule.
*   **Dependency Management**: All tool commands within the `shell` block of a rule should be prefixed with `pixi run` to ensure the use of the correct, version-locked tool.
