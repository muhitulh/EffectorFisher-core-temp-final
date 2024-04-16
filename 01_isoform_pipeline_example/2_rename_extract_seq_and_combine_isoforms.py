import csv
import glob

def get_isoform_name(index):
    return str(index + 1)

def extract_suffix(id_string):
    parts = id_string.split('_')
    return '_'.join(parts[1:3])+'_'

def clean_up_id(id_string):
    # Check if "_SNOO" is in the string
    if "_SNOO" in id_string:
        # Split the string at "_SNOO" and keep everything before it
        id_cleaned = id_string.split("_SNOO")[0]
        return id_cleaned
    # Return the original ID if "_SNOO" is not found
    return id_string

def process_file(input_file, output_file, mapping_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        mapping_writer = csv.writer(open(mapping_file, 'w', newline=''), delimiter='\t')
        
        original_header = next(reader, None)
        if not original_header:
            print(f"Warning: The file {input_file} is empty.")
            return

        first_data_row = next(reader, None)
        if not first_data_row:
            print(f"Warning: The file {input_file} has no data rows.")
            return
        suffix = extract_suffix(first_data_row[0])

        new_headers = [original_header[0]]
        mapping_headers = ['isoform', 'sequence']
        mapping_writer.writerow(mapping_headers)
        
        for i, col in enumerate(original_header[1:], start=1):
            new_name = suffix + get_isoform_name(i-1)
            new_headers.append(new_name)
            mapping_writer.writerow([new_name, col])
        
        writer.writerow(new_headers)
        
        # Clean up ID in the first_data_row before writing it
        first_data_row[0] = clean_up_id(first_data_row[0])
        writer.writerow(first_data_row)
        
        for row in reader:
            # Clean up ID for each row
            row[0] = clean_up_id(row[0])
            writer.writerow(row)

if __name__ == "__main__":
    input_directory = "01_input_from_metaeuk"
    input_files = glob.glob(f"{input_directory}/*_prot.fas.byseqid_table")
    
    for input_file in input_files:
        temp_output_file = input_file + '_renamed'
        mapping_file = input_file + '_mapping.csv'
        process_file(input_file, temp_output_file, mapping_file)
        print(f"Processed {input_file}, new file is {temp_output_file}, mapping file is {mapping_file}")



### combine gene specific isoform output



import glob
import pandas as pd
import os

# Define the input and output directories
input_directory = '01_input_from_metaeuk'
output_directory = '02_isoform_combined_table'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Find all files with the suffix "_prot.fas.byseqid_table_renamed" in the specified directory
file_paths = glob.glob(f'{input_directory}/*_prot.fas.byseqid_table_renamed')

# Initialize an empty DataFrame to store the merged data
merged_df = pd.DataFrame()

# Process each file and merge horizontally
for file_path in file_paths:
    # Read data from the current file
    df = pd.read_csv(file_path, delim_whitespace=True)

    # If the merged DataFrame is empty, set it as the current DataFrame
    if merged_df.empty:
        merged_df = df
    else:
        # Merge the data horizontally based on the ID column and sum the values
        merged_df = pd.merge(merged_df, df, on=df.columns[0], how='outer').fillna(0)

# After merging, ensure all numeric columns are integers
for column in merged_df.columns[1:]:  # Skip the first column assuming it's an identifier
    merged_df[column] = merged_df[column].astype(float).round().astype(int)

# Save the merged DataFrame to a new file in the output directory
output_file_path = os.path.join(output_directory, 'combined_isoform_output.txt')
merged_df.to_csv(output_file_path, sep='\t', index=False)
