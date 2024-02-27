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
output_file_path = os.path.join(output_directory, 'combined_isoform_output_global.txt')
merged_df.to_csv(output_file_path, sep='\t', index=False)
