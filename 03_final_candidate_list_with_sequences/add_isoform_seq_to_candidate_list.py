import pandas as pd

# Load the isoform sequences from 'isoform_seq.txt'
isoform_seq_df = pd.read_csv('isoform_seq.txt', sep='\t')
# Create a dictionary to map isoform names to sequences
isoform_seq_dict = dict(zip(isoform_seq_df['isoform'], isoform_seq_df['sequence']))

# Load the filtered locus list from '3_Filtered_locus_list_0.05.txt'
filtered_locus_df = pd.read_csv('candidate_list.txt', sep='\t') ##change input id here

# Map sequences to the locus list based on isoform names
filtered_locus_df['Sequence'] = filtered_locus_df['isoform'].map(isoform_seq_dict)

# Save the modified dataframe with sequences added to a new file
output_file_path = 'candidate_list_with_sequences.txt'
filtered_locus_df.to_csv(output_file_path, sep='\t', index=False)

print(f"Updated file saved to: {output_file_path}")
