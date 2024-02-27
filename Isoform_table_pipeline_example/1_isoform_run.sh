#!/bin/bash

# Change the path to the directory containing your files
input_directory="/mnt/c/Users/20435585/Desktop/Isoform_table_pipeline_example/01_input_from_metaeuk"

# Loop through all the files matching the pattern *_prot.fas in the input directory
for file in "$input_directory"/*_prot.fas; do
    echo "Running script for file: $file"
    perl cluster_exact_seqs_fasta.pl "$file"
done