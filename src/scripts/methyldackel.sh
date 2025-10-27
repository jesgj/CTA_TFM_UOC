#!/bin/bash

# --- User Configuration ---
# 1. Set the path to your reference genome FASTA file
GENOME_FASTA="$HOME/genomes/bluebottle/WGBS_phage/Physalia_utriculus.Dec23.genome.control_phage.fasta"

# 2. Set the path to the directory containing your BAM files
#    The script will look for a BAM file matching the sample name in this directory
BAM_DIR="results/filtered_bams/"

# 3. Set the path to the TSV file you just created
MBIAS_OPTIONS_FILE="results/mbias/all_samples_mbias_options.tsv"

# 4. Set the directory where you want the final methylKit output files to go
OUTPUT_DIR="results/methyldackel_out_methylKit/"
# --- End of Configuration ---


# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Read the TSV file, skipping the header line
tail -n +2 "${MBIAS_OPTIONS_FILE}" | while IFS=$'\t' read -r sample_name options; do
    # Construct the full path to the input BAM file
    # Note: This assumes the BAM file is named exactly as in the TSV, plus '.bam'
    input_bam="${BAM_DIR}/${sample_name}.bam"

    # Define the output prefix for this sample
    output_prefix="${OUTPUT_DIR}/${sample_name}"

    echo "--- Processing sample: ${sample_name} ---"

    # Check if the BAM file actually exists before trying to run the command
    if [ ! -f "${input_bam}" ]; then
        echo "ERROR: BAM file not found at ${input_bam}"
        echo "Skipping this sample."
        continue # Skip to the next iteration of the loop
    fi

    # Construct and execute the MethylDackel extract command.
    # The 'options' variable from the file is inserted directly into the command.
    # We use 'eval' to ensure that the shell correctly interprets the options string
    # with its multiple arguments.i
    # --methylkit not compatible with mergeContext (do it in a independent run)
    command="MethylDackel extract --methylKit "${options}" -o "${output_prefix}" --minOppositeDepth 10 --maxVariantFrac 0.5 "${GENOME_FASTA}" "${input_bam}""
    echo "Running command:"
    echo "${command}"
    eval "${command}"
    echo "Finished processing ${sample_name}."
    echo ""

done

echo "-------------------------------------------"
echo "All samples have been processed!"
echo "Your methylKit files are located in: ${OUTPUT_DIR}"
