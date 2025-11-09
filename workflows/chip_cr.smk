import os
import re
from collections import defaultdict

# --- CONFIGURATION ---
# Define directories from config
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
BOWTIE2_INDEX_DIR = config["bowtie2_index_dir"]
BAM_QC_DIR = config["bam_qc_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
DEEPTOOLS_DIR = config["deeptools_dir"]
BIGWIG_DIR = config["bigwig_dir"]
SUBTRACTED_BIGWIG_DIR = config["subtracted_bigwig_dir"]

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["alignment_dir"] = ALIGNMENT_DIR
config["bowtie2_index_dir"] = BOWTIE2_INDEX_DIR
config["bam_qc_dir"] = BAM_QC_DIR
config["filtered_bam_dir"] = FILTERED_BAM_DIR
config["filtered_bam_qc_dir"] = FILTERED_BAM_QC_DIR
config["deeptools_dir"] = DEEPTOOLS_DIR
config["bigwig_dir"] = BIGWIG_DIR
config["subtracted_bigwig_dir"] = SUBTRACTED_BIGWIG_DIR

# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(ALIGNMENT_DIR, exist_ok=True)
os.makedirs(BOWTIE2_INDEX_DIR, exist_ok=True)
os.makedirs(BAM_QC_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_QC_DIR, exist_ok=True)
os.makedirs(DEEPTOOLS_DIR, exist_ok=True)
os.makedirs(BIGWIG_DIR, exist_ok=True)
os.makedirs(SUBTRACTED_BIGWIG_DIR, exist_ok=True)
os.makedirs(os.path.join("logs", "chipseq_cutrun", "bowtie2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", "chipseq_cutrun", "bam_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", "chipseq_cutrun", "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join("logs", "chipseq_cutrun", "deeptools"), exist_ok=True)
os.makedirs(os.path.join("logs", "chipseq_cutrun", "bamCoverage"), exist_ok=True)


# --- SAMPLE DISCOVERY ---
def discover_samples_chip_cr(raw_dir):
    if not os.path.exists(raw_dir):
        return {}, []

    pe_files = defaultdict(dict)
    se_files = {}
    
    # First, find all potential PE files
    for f in os.listdir(raw_dir):
        if not f.endswith((".fq.gz", ".fastq.gz")):
            continue
        
        match = re.match(r"(.+?)_(?:R)?([12])\.(fastq|fq)\.gz$", f)
        if match:
            sample = match.group(1)
            read = match.group(2) # from "1" or "2"
            pe_files[sample][read] = os.path.join(raw_dir, f)
        else:
            # If it doesn't match PE pattern, it's a potential SE file
            match = re.match(r"(.+?)\.(fastq|fq)\.gz", f)
            if match:
                sample = match.group(1)
                se_files[sample] = os.path.join(raw_dir, f)

    samples_info = {}
    # Process PE files
    for sample, reads in pe_files.items():
        if '1' in reads and '2' in reads:
            samples_info[sample] = {'R1': reads['1'], 'R2': reads['2'], 'type': 'PE'}
            # If a file was mistakenly identified as SE, remove it
            if sample in se_files:
                del se_files[sample]

    # Process remaining as SE files
    for sample, path in se_files.items():
        # Make sure it's not a leftover from an incomplete PE pair
        if sample not in samples_info:
            samples_info[sample] = {'R1': path, 'R2': '', 'type': 'SE'}
            
    return samples_info, list(samples_info.keys())

SAMPLES_INFO, SAMPLES = discover_samples_chip_cr(RAW_DIR)

config['samples_info'] = SAMPLES_INFO

def get_subtraction_pairs(samples_info, bigwig_dir):
    pairs = []
    
    chip_samples = [s for s in samples_info if '_input' not in s]
    input_samples = [s for s in samples_info if '_input' in s]

    for chip_sample in chip_samples:
        # Construct the expected input name, e.g., "liver_rep1" -> "liver_input_rep1"
        parts = chip_sample.split('_')
        if len(parts) == 2: # Expects "tissue_repX" format
            base, rep = parts
            expected_input = f"{base}_input_{rep}"
            
            if expected_input in input_samples:
                chip_info = samples_info[chip_sample]
                input_info = samples_info[expected_input]
                
                chip_read_type = 'pe' if chip_info['type'] == 'PE' else 'se'
                # Input and chip might have different read types, but for pairing we assume they match
                # The bw file will have the correct read_type suffix from its own rule
                input_read_type = 'pe' if input_info['type'] == 'PE' else 'se'

                chip_bw = os.path.join(bigwig_dir, f"{chip_sample}_{chip_read_type}.bw")
                input_bw = os.path.join(bigwig_dir, f"{expected_input}_{input_read_type}.bw")
                
                pairs.append({
                    'chip_bw': chip_bw, 
                    'input_bw': input_bw, 
                    'sample': chip_sample, 
                    'read_type': chip_read_type
                })
                
    return pairs

SUBTRACTION_PAIRS = get_subtraction_pairs(SAMPLES_INFO, BIGWIG_DIR)

# Pass SUBTRACTION_PAIRS to config so it's accessible in included rules
config['subtraction_pairs'] = SUBTRACTION_PAIRS

# --- MODULE INCLUSION ---

# 1. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"

# 3. Alignment
include: "rules/chip_cr/alignment.smk"

# 4. BAM QC
include: "rules/chip_cr/bam_qc.smk"

# 5. Filtering and Deduplication
include: "rules/chip_cr/filter_bam.smk"

# 6. Filtered BAM QC (QC + deeptools summary)
include: "rules/chip_cr/filtered_bam_qc.smk"

# 7. BigWig generation and subtraction
include: "rules/chip_cr/bigwig.smk"

# 8. Heatmap computeMatrix + plotHeatmap
include: "rules/chip_cr/heatmap.smk"


# --- FINAL TARGETS ---

def get_all_final_outputs(samples_info):
    outputs = []

    # Add raw QC outputs
    for sample, info in samples_info.items():
        if info['type'] == 'PE':
            outputs.append(os.path.join(QC_DIR, f"{sample}_R1_raw_fastqc.html"))
            outputs.append(os.path.join(QC_DIR, f"{sample}_R2_raw_fastqc.html"))
        else: # SE
            outputs.append(os.path.join(QC_DIR, f"{sample}_raw_fastqc.html"))

    # Add trimmed QC outputs
    for sample, info in samples_info.items():
        if info['type'] == 'PE':
            outputs.append(os.path.join(QC_TRIMMED_DIR, f"{sample}_R1_trimmed_fastqc.html"))
            outputs.append(os.path.join(QC_TRIMMED_DIR, f"{sample}_R2_trimmed_fastqc.html"))
        else: # SE
            outputs.append(os.path.join(QC_TRIMMED_DIR, f"{sample}_SE_trimmed_fastqc.html"))

    # Add BAM QC outputs
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        outputs.extend([
            os.path.join(BAM_QC_DIR, f"{sample}_{read_type}.stats.txt"),
            os.path.join(BAM_QC_DIR, f"{sample}_{read_type}.flagstat.txt"),
            os.path.join(BAM_QC_DIR, f"{sample}_{read_type}.alignment_summary_metrics.txt")
        ])

    # Add filtered BAM files
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        outputs.append(os.path.join(FILTERED_BAM_DIR, f"{sample}_{read_type}.filtered.sorted.bam"))

    # Add filtered BAM QC outputs
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        outputs.extend([
            os.path.join(FILTERED_BAM_QC_DIR, f"{sample}_{read_type}.stats.txt"),
            os.path.join(FILTERED_BAM_QC_DIR, f"{sample}_{read_type}.flagstat.txt"),
            os.path.join(FILTERED_BAM_QC_DIR, f"{sample}_{read_type}.alignment_summary_metrics.txt")
        ])

    # Add deeptools outputs
    outputs.extend([
        os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"),
        os.path.join(DEEPTOOLS_DIR, "fingerprints.png"),
        os.path.join(DEEPTOOLS_DIR, "heatmap.png")
    ])

    # Add bigwig outputs
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        outputs.append(os.path.join(BIGWIG_DIR, f"{sample}_{read_type}.bw"))

    # Add subtracted bigwig outputs
    for pair in config['subtraction_pairs']:
        outputs.append(os.path.join(SUBTRACTED_BIGWIG_DIR, f"{pair['sample']}_{pair['read_type']}.subtracted.bw"))

    return outputs

rule all:
    input:
        get_all_final_outputs(SAMPLES_INFO)
