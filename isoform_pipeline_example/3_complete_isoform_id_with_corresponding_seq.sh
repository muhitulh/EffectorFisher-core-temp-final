#!/bin/bash

# Specify input and output directories
input_directory="01_input_from_metaeuk"
output_directory="03_isoform_seq"

# Check if the output directory exists; if not, create it
if [ ! -d "$output_directory" ]; then
    mkdir -p "$output_directory"
fi

# Construct the output file path
output_file="${output_directory}/isoform_seq.txt"

# Get the first file from the input directory and include its headers in the output file
head -n 1 "$(ls ${input_directory}/*.csv | head -n 1)" > "$output_file"

# Concatenate all files from the input directory but skip the headers
for f in ${input_directory}/*.csv; do
    tail -n +2 "$f" >> "$output_file"
done
