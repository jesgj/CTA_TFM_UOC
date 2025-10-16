# workflows/rules/bam_qc.smk
import os

# --- CONFIGURATION ---
DEDUP_DIR = config["dedup_dir"]
DEDUP_BAM_QC_DIR = config["dedup_bam_qc_dir"]
REF_GENOME = config["ref_genome"]
SAMPLES = list(config['samples_info'].keys())

# --- RULES ---

rule samtools_stats_dedup:
    """
    Generate alignment statistics for a deduplicated BAM file using samtools stats.
    """
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        stats=os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.stats.txt")
    threads: 4
    log:
        os.path.join("logs", "samtools_stats_dedup", "{sample}.log")
    shell:
        """
        pixi run samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}
        """

rule samtools_flagstat_dedup:
    """
    Generate flag statistics for a deduplicated BAM file using samtools flagstat.
    """
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        flagstat=os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.flagstat.txt")
    threads: 1
    log:
        os.path.join("logs", "samtools_flagstat_dedup", "{sample}.log")
    shell:
        """
        pixi run samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        """

rule picard_collect_alignment_metrics_dedup:
    """
    Collects alignment metrics from a deduplicated BAM file using Picard Tools.
    """
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam"),
        ref=REF_GENOME
    output:
        metrics=os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.alignment_summary_metrics.txt")
    params:
        # Picard needs a lot of memory
        java_opts = "-Xmx4g"
    log:
        os.path.join("logs", "picard_collect_alignment_metrics_dedup", "{sample}.log")
    shell:
        """
        pixi run picard {params.java_opts} CollectAlignmentSummaryMetrics \
            R={input.ref} \
            I={input.bam} \
            O={output.metrics} > {log} 2>&1
        """

# Rule to aggregate all deduplicated BAM QC results
rule all_dedup_bam_qc:
    input:
        expand(os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.stats.txt"), sample=SAMPLES),
        expand(os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.flagstat.txt"), sample=SAMPLES),
        expand(os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.alignment_summary_metrics.txt"), sample=SAMPLES)
