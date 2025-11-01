# workflows/rules/chip_cr/filtered_bam_qc.smk
import os

# --- CONFIGURATION ---
FILTERED_BAM_DIR = config["filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
REF_GENOME = config["ref_genome"]
DEEPTOOLS_DIR = config["deeptools_dir"]
SAMPLES_INFO = config['samples_info']
BIGWIG_DIR = config["bigwig_dir"]
SUBTRACTED_BIGWIG_DIR = config["subtracted_bigwig_dir"]

# Get extra args from config (Deeptools)
MBS_ARGS = config.get("deeptools", {}).get("multiBamSummary", {}).get("extra_args", "")
PC_ARGS = config.get("deeptools", {}).get("plotCorrelation", {}).get("extra_args", "")
PF_ARGS = config.get("deeptools", {}).get("plotFingerprint", {}).get("extra_args", "")

# --- HELPER FUNCTION ---
def get_filtered_bam(wildcards):
    """
    Returns the path to the filtered BAM file for a sample,
    using the {read_type} wildcard ('pe' or 'se').
    """
    return os.path.join(FILTERED_BAM_DIR, f"{wildcards.sample}_{wildcards.read_type}.filtered.sorted.bam")

def get_all_filtered_bams(samples_info):
    bams = []
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        bams.append(os.path.join(FILTERED_BAM_DIR, f"{sample}_{read_type}.filtered.sorted.bam"))
    return bams

ALL_FILTERED_BAMS = get_all_filtered_bams(SAMPLES_INFO)

# --- QC RULES (Samtools & Picard) ---

rule samtools_stats_filtered:
    """
    Generate alignment statistics for a filtered BAM file using samtools stats.
    """
    input:
        bam = get_filtered_bam
    output:
        stats = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_{read_type}.stats.txt")
    threads: 4
    log:
        os.path.join("logs", "chipseq_cutrun", "bam_qc", "{sample}_{read_type}_filtered_stats.log")
    shell:
        "pixi run samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}"

rule samtools_flagstat_filtered:
    """
    Generate flag statistics for a filtered BAM file using samtools flagstat.
    """
    input:
        bam = get_filtered_bam
    output:
        flagstat = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_{read_type}.flagstat.txt")
    threads: 1
    log:
        os.path.join("logs", "chipseq_cutrun", "bam_qc", "{sample}_{read_type}_filtered_flagstat.log")
    shell:
        "pixi run samtools flagstat {input.bam} > {output.flagstat} 2> {log}"

rule picard_collect_alignment_metrics_filtered:
    """
    Collects alignment metrics from a filtered BAM file using Picard Tools.
    """
    input:
        bam = get_filtered_bam,
        ref = REF_GENOME
    output:
        metrics = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_{read_type}.alignment_summary_metrics.txt")
    params:
        java_opts = "-Xmx4g"
    log:
        os.path.join("logs", "chipseq_cutrun", "bam_qc", "{sample}_{read_type}_filtered_picard_metrics.log")
    shell:
        """
pixi run picard {params.java_opts} CollectAlignmentSummaryMetrics \
            R={input.ref} \
            I={input.bam} \
            O={output.metrics} > {log} 2>&1
        """

# --- DEEPTOOLS RULES ---

rule multiBamSummary:
    """
    Computes read coverages for multiple BAM files.
    """
    input:
        bams = ALL_FILTERED_BAMS
    output:
        npz = os.path.join(DEEPTOOLS_DIR, "read_coverage.npz")
    params:
        extra = MBS_ARGS
    threads: 8
    log:
        os.path.join("logs", "chipseq_cutrun", "deeptools", "multiBamSummary.log")
    shell:
        "pixi run multiBamSummary bins -b {input.bams} -o {output.npz} -p {threads} {params.extra} &> {log}"

rule plotCorrelation:
    """
    Creates a heatmap of correlations based on multiBamSummary output.
    """
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "read_coverage.npz")
    output:
        heatmap = os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"),
        matrix = os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab")
    params:
        extra = PC_ARGS
    threads: 1
    log:
        os.path.join("logs", "chipseq_cutrun", "deeptools", "plotCorrelation.log")
    shell:
        "pixi run plotCorrelation -in {input.npz} -o {output.heatmap} --outFileCorMatrix {output.matrix} {params.extra} &> {log}"

rule plotFingerprint:
    """
    Generates fingerprints for each BAM file to assess ChIP-seq quality.
    """
    input:
        bams = ALL_FILTERED_BAMS
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "fingerprints.png"),
        metrics = os.path.join(DEEPTOOLS_DIR, "fingerprints.metrics.tab")
    params:
        extra = PF_ARGS
    threads: 8
    log:
        os.path.join("logs", "chipseq_cutrun", "deeptools", "plotFingerprint.log")
    shell:
        "pixi run plotFingerprint -b {input.bams} -o {output.plot} --outRawCounts {output.metrics} -p {threads} {params.extra} &> {log}"

# (BigWig generation and heatmap rules were extracted to
#  'workflows/rules/chip_cr/bigwig.smk' and 'workflows/rules/chip_cr/heatmap.smk')
