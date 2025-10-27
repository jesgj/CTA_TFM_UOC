# workflows/rules/filtered_bam_qc.smk
import os

# --- CONFIGURATION ---
SORTED_FILTERED_BAM_DIR = config["sorted_filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
REF_GENOME = config["ref_genome"]
SAMPLES = list(config['samples_info'].keys())

# --- RULES ---

rule samtools_stats_filtered:
    """
    Generate alignment statistics for a filtered BAM file using samtools stats.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        stats=os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.stats.txt")
    threads: 4
    log:
        os.path.join("logs", "samtools_stats_filtered", "{sample}.log")
    shell:
        """
        pixi run samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}
        """

rule samtools_flagstat_filtered:
    """
    Generate flag statistics for a filtered BAM file using samtools flagstat.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        flagstat=os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.flagstat.txt")
    threads: 1
    log:
        os.path.join("logs", "samtools_flagstat_filtered", "{sample}.log")
    shell:
        """
        pixi run samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        """

rule picard_collect_alignment_metrics_filtered:
    """
    Collects alignment metrics from a filtered BAM file using Picard Tools.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
        ref=REF_GENOME
    output:
        metrics=os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.alignment_summary_metrics.txt")
    params:
        # Picard needs a lot of memory
        java_opts = "-Xmx4g"
    log:
        os.path.join("logs", "picard_collect_alignment_metrics_filtered", "{sample}.log")
    shell:
        """
        pixi run picard {params.java_opts} CollectAlignmentSummaryMetrics \
            R={input.ref} \
            I={input.bam} \
            O={output.metrics} > {log} 2>&1
        """

# Rule to aggregate all filtered BAM QC results
rule all_filtered_bam_qc:
    input:
        expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.stats.txt"), sample=SAMPLES),
        expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.flagstat.txt"), sample=SAMPLES),
        expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.alignment_summary_metrics.txt"), sample=SAMPLES)
