# workflows/rules/chip_cr/bam_qc.smk
import os

# --- CONFIGURATION ---
ALIGNMENT_DIR = config["alignment_dir"]
BAM_QC_DIR = config["bam_qc_dir"]
REF_GENOME = config["ref_genome"]

# --- HELPER FUNCTION ---
def get_aligned_bam(wildcards):
    """
    Returns the path to the aligned BAM file for a sample,
    using the {read_type} wildcard ('pe' or 'se').
    """
    return os.path.join(ALIGNMENT_DIR, f"{wildcards.sample}_{wildcards.read_type}.sorted.bam")

# --- QC RULES ---

rule samtools_stats:
    """
    Generate alignment statistics for a BAM file using samtools stats.
    """
    input:
        bam = get_aligned_bam
    output:
        stats = os.path.join(BAM_QC_DIR, "{sample}_{read_type}.stats.txt")
    threads: 4
    log:
        os.path.join("logs", "chipseq_cutrun", "bam_qc", "{sample}_{read_type}_stats.log")
    shell:
        "pixi run samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}"

rule samtools_flagstat:
    """
    Generate flag statistics for a BAM file using samtools flagstat.
    """
    input:
        bam = get_aligned_bam
    output:
        flagstat = os.path.join(BAM_QC_DIR, "{sample}_{read_type}.flagstat.txt")
    threads: 1
    log:
        os.path.join("logs", "chipseq_cutrun", "bam_qc", "{sample}_{read_type}_flagstat.log")
    shell:
        "pixi run samtools flagstat {input.bam} > {output.flagstat} 2> {log}"

rule picard_collect_alignment_metrics:
    """
    Collects alignment metrics from a BAM file using Picard Tools.
    """
    input:
        bam = get_aligned_bam,
        ref = REF_GENOME
    output:
        metrics = os.path.join(BAM_QC_DIR, "{sample}_{read_type}.alignment_summary_metrics.txt")
    params:
        java_opts = "-Xmx4g"
    log:
        os.path.join("logs", "chipseq_cutrun", "bam_qc", "{sample}_{read_type}_picard_metrics.log")
    shell:
        """
        pixi run picard {params.java_opts} CollectAlignmentSummaryMetrics \
            R={input.ref} \
            I={input.bam} \
            O={output.metrics} > {log} 2>&1
        """
