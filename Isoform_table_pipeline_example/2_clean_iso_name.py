import glob
import os

# Define the pattern to match all input files
input_files_pattern = os.path.join('01_input_from_metaeuk', '*_prot.fas.byseqid_table')
# Define the base directory for the output files, to save in the current folder or another specific folder
output_base_dir = '01_input_from_metaeuk'  # Adjust this to your needs

# Use glob to find all files matching the input_files_pattern
input_files = glob.glob(input_files_pattern)

for original_file_path in input_files:
    # Construct the final output file path based on the original file name, adding '_renamed' before the file extension
    base_name = os.path.basename(original_file_path)
    output_file_name = base_name.replace('_prot.fas.byseqid_table', '_prot.fas.byseqid_table_renamed')
    output_file_path = os.path.join(output_base_dir, output_file_name)
    
    # Define the path for the file to save headers comparison, also based on the original file name
    headers_comparison_file_name = base_name.replace('_prot.fas.byseqid_table', '_isoform_seq.txt')
    headers_comparison_file_path = os.path.join(output_base_dir, headers_comparison_file_name)

    # Extract the first three parts of the original file name for column naming
    file_name_parts = base_name.split('_')[:3]
    col_name_base = '_'.join(file_name_parts)

    # Ensure the output directory exists
    os.makedirs(output_base_dir, exist_ok=True)

    # Process each file
    with open(original_file_path, 'r') as infile, open(output_file_path, 'w') as outfile, open(headers_comparison_file_path, 'w') as headers_file:
        for line_num, line in enumerate(infile):
            parts = line.rstrip().split("\t")
            # Modify the ID in the first column to keep only the part before the first underscore
            parts[0] = parts[0].split('_')[0]
            
            if line_num == 0:  # If this is the header row
                # Store old headers
                old_headers = parts
                # Modify column names starting from the second column
                new_headers = [parts[0]] + [f"{col_name_base}_{i}" for i in range(1, len(parts))]
                outfile.write("\t".join(new_headers) + "\n")
                
                # Write headers comparison
                headers_file.write("Old Header\tNew Header\n")
                for old, new in zip(old_headers, new_headers):
                    headers_file.write(f"{old}\t{new}\n")
            else:
                outfile.write("\t".join(parts) + "\n")