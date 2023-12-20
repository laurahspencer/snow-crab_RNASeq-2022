#!/bin/bash

# Input and output file paths
input_file="input.fasta"
output_file="output.fasta"

# Number of contigs to merge in each group
contigs_per_group=100

# Number of N's to separate merged sequences
separator_length=1000

# Temporary file for storing intermediate results
temp_file="temp.fasta"

# Create temporary file
touch "$temp_file"

# Count the total number of contigs in the input file
total_contigs=$(grep -c "^>" "$input_file")

# Calculate the number of groups
num_groups=$((total_contigs / contigs_per_group))

# Loop through each group
for ((i = 1; i <= num_groups; i++)); do
    # Calculate start and end indices for each group
    start_index=$(( (i - 1) * contigs_per_group + 1 ))
    end_index=$(( i * contigs_per_group ))

    # Extract contigs for the current group
    sed -n "${start_index},${end_index}p" "$input_file" > "$temp_file"

    # Merge contigs with N's separator
    merged_sequence=$(awk 'BEGIN {ORS="N"} NR>1 {print ""} {printf "%s", $0} END {print ""}' "$temp_file")

    # Append merged sequence to the output file
    echo -e ">$i\n$merged_sequence" >> "$output_file"
done

# Clean up temporary files
rm "$temp_file"

echo "Merging complete. Output written to $output_file"
