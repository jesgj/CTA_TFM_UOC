#!/bin/bash

# Set the path to your reference genome FASTA file
GENOME_FASTA="$HOME/genomes/bluebottle/WGBS_phage/Physalia_utriculus.Dec23.genome.control_phage.fasta"

# Set the path to the directory containing your BAM files
BAM_DIR="results/filtered_bams/"

OUTPUT_DIR="results/mbias/"


# Set the name for the final output TSV file
OUTPUT_FILE="${OUTPUT_DIR}"/"all_samples_mbias_options.tsv"
# --- End of Configuration ---

mkdir -p $OUTPUT_DIR

# Initialize the output file with a header
# The '-e' flag enables interpretation of backslash escapes like '\t' for tab
echo -e "sample\toptions" > "${OUTPUT_FILE}"

# Loop through every file ending with .bam in the specified directory
for bam_file in ${BAM_DIR}/*.bam; do
    # Get a clean sample name by removing the path and the .bam extension
    sample_name=$(basename "${bam_file}" .bam)

    echo "Processing sample: ${sample_name}"

    # Run MethylDackel mbias.
    # '2>&1' redirects stderr (where the message is) to stdout.
    # We then pipe '|' this output to 'grep' to find only the line we want.
    # 'sed' is then used to remove the descriptive part of the line.
    # The final result (just the options) is stored in the 'options' variable.
    options=$(MethylDackel mbias -@ 16 "${GENOME_FASTA}" "${bam_file}" "${OUTPUT_DIR}/${sample_name}" 2>&1 | \
              grep "Suggested inclusion options:" | \
              sed 's/Suggested inclusion options: //')

    # Check if we successfully captured the options
    if [ -n "${options}" ]; then
        # Append the sample name and the captured options, separated by a tab, to our output file
        echo -e "${sample_name}\t${options}" >> "${OUTPUT_FILE}"
    else
        echo "WARNING: Could not find suggested options for ${sample_name}."
    fi
done

echo "-------------------------------------------"
echo "Script finished!"
echo "Your TSV file is ready: ${OUTPUT_FILE}"
echo ""
echo "--- Contents of ${OUTPUT_FILE} ---"
cat "${OUTPUT_FILE}"
