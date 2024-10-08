#!/usr/bin/env bash
set -e
set -u

# Specify the path for the input files and the output files
reference_path=$1
input_path=$2
output_path=$3

# Set the number of threads to use
threads=$4

# Read the reference sequence
reference=$(awk 'NR==2' $reference_path)

# Ensure the output path exists
mkdir -p $output_path

# Loop through all the reads fastq files in the input path
for reads_file in $input_path/reads_fc*.fastq.gz; do
    # Extract the number after fc
    number=$(basename $reads_file | sed -E 's/reads_fc([0-9]+)\.fastq.gz/\1/')

    # Generate uuid
    uuid=$(uuidgen)
  
    # Generate the corresponding seed file and output bam file names
    reference_file="$input_path/reference_fc${number}.fasta"
    output_bam="$output_path/calls_to_draft_fc${number}"

    echo "Processing reads file: $reads_file."
    mini_align.sh -i $reads_file -r $reference_file -m -p $output_bam -t $threads

done