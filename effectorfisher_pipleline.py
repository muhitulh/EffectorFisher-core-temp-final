## for thresholds set up in command line
#usage: python 0_effectorfisher_pipleline.py --min_iso 5 --cyst 2 --pred_score 2 --total_aa 300 --p_value 0.05


import pandas as pd
import argparse
import glob
import os
import math
import csv

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process isoform frequencies.')
# Add this line in the section where you setup argparse.ArgumentParser
parser.add_argument('--min_iso', type=int, default=10, help='Minimum isoform count threshold')
parser.add_argument('--cyst', type=float, default=2, help='Cysteine count threshold')
parser.add_argument('--pred_score', type=float, default=2, help='Prediction score threshold')
parser.add_argument('--total_aa', type=float, default=300, help='Total amino acid count threshold')
parser.add_argument('--p_value', type=float, default=0.05, help='P-value threshold')

# Parse arguments
args = parser.parse_args()

###Step 3: calculate median and assign disease to low if value lower than median, and high - if value higher than median---------------------------
#import pandas as pd
#import glob

# Using glob to find all files starting with '1_data' in the specified directory
data_files = glob.glob('1_data*.txt')

# Processing each file
for data_file in data_files:
    data_df = pd.read_csv(data_file, sep='\t')

    # Check if 'disease' column exists

    
    if 'disease' in data_df.columns:
        # Calculate the median
        median_value = data_df['disease'].median()

        # Replace values based on the median
        data_df['disease'] = data_df['disease'].apply(lambda x: 'low' if x < median_value else 'high')

        # Save the modified dataframe back to the file
        data_df.to_csv(data_file, sep='\t', index=False)
    else:
        print(f"'disease' column not found in {data_file}")
        
print("step 3 completed - set disease severity level.")
        

## Step 1 (&2): create cultivar specific files--------------------------original step 2 deleted, as that included here------------------------

#import pandas as pd

# Load the file, ensuring that the first column is treated as the 'ID' column
df = pd.read_csv('input_files/0_phenotype_data.txt', sep='\t')

# If the dataframe does not have headers, you need to add them
if df.columns[0] != 'ID':
    # Assuming the first column is the one to be renamed to 'ID'
    df.columns = ['ID'] + df.columns[1:].tolist()

# Splitting the file and changing the second column header to 'disease'
for i, col_name in enumerate(df.columns[1:], start=1):  # start=1 to save files with correct cultivar number
    # Create a subset dataframe with just the 'ID' column and the current cultivar's data
    subset_df = df[['ID', col_name]].copy()
    
    # Rename the second column to 'disease'
    subset_df.columns = ['ID', 'disease']
    
    # Save this dataframe to a file with the original column name in the file name
    output_file_name = f'1_data{i}.txt'
    subset_df.to_csv(output_file_name, sep='\t', index=False)

print("step 1 completed: cultivar-specific files with 'disease' column have been created.")

## Step 2: remove rows (isolates) with missing disease data-----------------------

#import pandas as pd
#import glob

# Using glob to find all files starting with '1_data' in the specified directory
data_files = glob.glob('1_data*.txt')

# Process each file
for data_file in data_files:
    # Read the data file
    data_df = pd.read_csv(data_file, sep='\t')

    # Remove rows where the second column is empty
    # Assuming the second column is at index 1
    data_df = data_df[data_df.iloc[:, 1].notna()]

    # Optionally, save the modified DataFrame back to the file
    data_df.to_csv(data_file, sep='\t', index=False)

    # Print the name of the processed file (optional)
    print(f"Processed {data_file}")

print("step 2 completed: rows with missing disease data removed.")




## Step 4: process the complete isoform table - remove isoform frequency <5--------------------------------------------
##usage: python step_4_effectorfisher_pipeline.py --min-iso 10 
# --min-iso default value is 5
#import pandas as pd
#import argparse

# Set up the argument parser
#parser = argparse.ArgumentParser(description='Process isoform frequencies.')
#parser.add_argument('--min-iso', type=int, default=5, help='Minimum isoform frequency for inclusion')

# Parse the command-line arguments
#args = parser.parse_args()

# Load the file into a pandas DataFrame
df = pd.read_csv('input_files/0_combined_isoform.txt', sep='\t')

# Sum the frequencies of each isoform across all samples
isoform_sums = df.iloc[:, 1:].sum()

# Filter out isoforms with a sum less than the --min-iso value
filtered_isoforms = isoform_sums[isoform_sums > args.min_iso]

# Creating a new DataFrame with only the filtered isoforms
filtered_df = df[['ID'] + list(filtered_isoforms.index)]

# Saving the filtered DataFrame to a new file
filtered_df.to_csv('input_files/0_filtered_combined_isoform.txt', index=False, sep='\t')

print("step 4 completed")

## Step 5: change 1.0 = P and 0.0 = A----------------------------------------------------------------------------
import pandas as pd
import glob

# Using glob to find all files starting with '1_data' in the specified directory
data_files = glob.glob('input_files/0_filtered_combined_isoform.txt')

# Processing each file
for data_file in data_files:
    data_df = pd.read_csv(data_file, sep='\t')

    # Looping through all columns after the second one
    for col in data_df.columns[1:]:
        # Replace 1 with 'P' and 0 with 'A'
        data_df[col] = data_df[col].apply(lambda x: 'P' if x == 1 else ('A' if x == 0 else x))

    # Save the modified dataframe back to the file
    data_df.to_csv(data_file, sep='\t', index=False)
print("step 5 completed")



####step 6: merge cultivar-specific file with combine isoform file by ID--------------------------------------------
#import pandas as pd
#import glob
#import os

# File to be concatenated
additional_file = 'input_files/0_filtered_combined_isoform.txt'

# Read the additional file
additional_df = pd.read_csv(additional_file, sep='\t')

# Using glob to find all files starting with '1_data' in the specified directory
data_files = glob.glob( '1_data*.txt')

# Concatenating each file with the additional file
for data_file in data_files:
    # Read the data file
    data_df = pd.read_csv(data_file, sep='\t')

    # Concatenate with the additional dataframe on 'ID' column
    concatenated_df = pd.merge(data_df, additional_df, on='ID', how='left')

    # Corrected extraction of file number
    file_number = os.path.basename(data_file).replace('1_data', '').replace('.txt', '')
    
    # Save the concatenated dataframe to a new file
    concatenated_file_path =  f'2_merged_data{file_number}.txt'
    concatenated_df.to_csv(concatenated_file_path, sep='\t', index=False)
    
print("step 6 completed")
    
    
## step 7: contingency table------------------------------------------
#import pandas as pd
#import glob

# Find all files that start with '4_merged_' and end with '.txt'
file_names = glob.glob('2_merged_*.txt')

for file_name in file_names:
    # Load the data
    data = pd.read_csv(file_name, delimiter='\t')

    # Automatically identify isoform columns (all columns starting from the third one)
    isoform_columns = data.columns[2:]

    # Initialize a dictionary to hold the counts for each category
    counts = {('high', 'A'): {}, ('high', 'P'): {}, ('low', 'A'): {}, ('low', 'P'): {}}

    # Initialize counts for each isoform in each category
    for isoform in isoform_columns:
        for key in counts.keys():
            counts[key][isoform] = 0

    # Iterate over each row to update counts
    for _, row in data.iterrows():
        disease_status = row['disease'].split()[0]  # Extract 'high' or 'low' from 'disease' column
        for isoform in isoform_columns:
            value = row[isoform]  # 'A' or 'P'
            if pd.notna(value):  # Make sure to handle NaN values
                counts[(disease_status, value)][isoform] += 1

    # Convert the counts dictionary to a DataFrame for easy handling
    contingency_df = pd.DataFrame.from_dict({
        (disease, value): counts[(disease, value)]
        for disease in ['high', 'low']
        for value in ['A', 'P']
    }, orient='index').transpose()

    # Rename the columns to match the desired output format
    contingency_df.columns = ['{}-{}'.format(*col) for col in contingency_df.columns]

    # Reorder the columns according to the specified order if necessary
    desired_column_order = sorted(contingency_df.columns, key=lambda x: (x.split('-')[1], x.split('-')[0]))
    contingency_df = contingency_df[desired_column_order]

    # Define the output file name based on the input file name
    output_file_name = file_name.replace('2_merged_', '3_contingency_')
    # Save the corrected contingency table to a file
    contingency_df.to_csv(output_file_name, sep='\t', index_label='isoform')

    print(f"Processed and saved contingency table: {output_file_name}")

    # Rename the columns for hypergeometric analysis
    column_rename = {'high-A': 'c', 'high-P': 'a', 'low-A': 'd', 'low-P': 'b'}
    contingency_df.rename(columns=column_rename, inplace=True)

    # Write the reformatted contingency table to a new file for hypergeometric analysis
    hypergeo_data_file = file_name.replace('2_merged_', '4_hypergeo_')
    contingency_df.to_csv(hypergeo_data_file, sep='\t', index=True, index_label='isoform')

    print(f"Processed and saved hypergeometric analysis table: {hypergeo_data_file}")

print("step 7 completed: contingency table")
    
    
## Step 8: hypergeomatric test-------------------------------------------------
###Hypergeo test - loop over many datasets together [automated]
###loop over all dataset
#import math
#import glob

# Find all files starting with '6_hypergeo_'
input_files = glob.glob('4_hypergeo_*.txt')

# Factorial calculation function
def calculate_factorials(n):
    fact = [0] * (n + 1)
    fact[0] = 1
    for i in range(2, n + 1):
        fact[i] = fact[i - 1] + math.log(i)
    return fact

def process_dataset(input_file, output_file, fact):
    with open(output_file, 'w') as fileout, open(input_file, 'r') as input_data:
        fileout.write("isoform\thigh-A\thigh-P\tlow-A\tlow-P\tp-value\n")
        next(input_data)  # Skip the header row
        for line_num, line in enumerate(input_data, 2):
            tabdata = line.strip().split('\t')
            try:
                iso = tabdata[0]
                a = round(float(tabdata[1]))
                c = round(float(tabdata[2]))
                b = round(float(tabdata[3]))
                d = round(float(tabdata[4]))
            except ValueError as e:
                print(f"Error on line {line_num}: {e}")
                continue
            
            n = a + b + c + d
            hyp1 = fact[a + b] + fact[c + d] + fact[a + c] + fact[b + d] - (fact[n] + fact[a] + fact[b] + fact[c] + fact[d])
            hyp1 = math.exp(hyp1)
            fileout.write(f"{iso}\t{a}\t{c}\t{b}\t{d}\t{hyp1}\n")

# Assuming calculate_factorials and process_dataset are predefined functions
# Pre-calculate factorials
fact = calculate_factorials(1000000)

# Loop over the input files
for input_file in input_files:
    # Extract the dataset number or name
    base_name = input_file.replace('4_hypergeo_', '')  # Remove '6_hypergeo_' prefix
    
    # Construct the output file name
    output_file = f'5_output_hypergeo_{base_name}'
    
    print(f"Processing {input_file}...")
    process_dataset(input_file, output_file, fact)

print("step 8 completed: hypergeomatric test")


#step 9 (&10) merge all the cultivar-specific p-value (master data) and add lowest p-value col-----------------------------------------------------
#import pandas as pd
#import glob

# Specify the pattern for your files
file_pattern = '5_output_hypergeo_data*.txt'
files = glob.glob(file_pattern)

# Initialize an empty DataFrame for the merged data
merged_data = pd.DataFrame()

# Process each file
for file_number, file_path in enumerate(sorted(files), start=1):
    # Read the file, focusing only on the 'isoform' and 'p-value' columns
    data = pd.read_csv(file_path, delimiter='\t', usecols=['isoform', 'p-value'])
    
    # Rename the 'p-value' column to 'p-value-{i}'
    p_value_column = f'p-value-{file_number}'
    data.rename(columns={'p-value': p_value_column}, inplace=True)
    
    # If merged_data is empty, initialize it with the data from the first file
    if merged_data.empty:
        merged_data = data
    else:
        # Merge the current data with the merged_data on 'isoform', using an outer join
        merged_data = pd.merge(merged_data, data, on='isoform', how='outer')

# After merging, calculate the lowest p-value across the merged 'p-value' columns
p_value_cols = [col for col in merged_data.columns if col.startswith('p-value')]
merged_data["p-value_lowest"] = merged_data[p_value_cols].min(axis=1)

# Save the merged data with the lowest p-value column to a new file
merged_data.to_csv('6_merged_lowest_p-value.txt', sep='\t', index=False)

print("Step 9 completed. Merged data with lowest p-value written to 'merged_lowest_p-value.txt'")



##step 10: need to create col locus_id from isoform_id---------------------------------------
###preparing cultivar specific fisher p-value for combining
#import pandas as pd

# Read the merged output
merged_df = pd.read_csv("6_merged_lowest_p-value.txt", delimiter="\t")

# Duplicate "isoform_id" column for further processing
merged_df["locus_id"] = merged_df["isoform"]

# Use regex to process the "variable" column
merged_df["locus_id"] = merged_df["locus_id"].str.replace(r'(_\d+)$', '', regex=True)

# Use case-insensitive regex to remove A, B, C, D, or any combination thereof from the end
merged_df["locus_id"] = merged_df["locus_id"].str.replace(r'(?i)([A-D]+)$', '', regex=True)

# Move "variable" column to the first position
cols = ['locus_id'] + [col for col in merged_df if col != 'locus_id']
merged_df = merged_df[cols]

# Save the updated dataframe
merged_df.to_csv("7_merged_lowest_p-value_with_locus_id.txt", index=False, sep="\t")

print("Step 10 completed: locus id added'")



###step 11: combine fisher result with predector result----------------------------------------(try with raw predector output)
#import pandas as pd

# Load the data files
df1 = pd.read_csv("7_merged_lowest_p-value_with_locus_id.txt", sep="\t")
df2 = pd.read_csv("input_files/0_predector_results.txt", sep="\t")

# Ensure that 'locus_id' is the first column in both dataframes
df1 = df1[['locus_id'] + [col for col in df1.columns if col != 'locus_id']]
df2 = df2[['locus_id'] + [col for col in df2.columns if col != 'locus_id']]

# Merge the dataframes on 'locus_id'
merged_df = df1.merge(df2, on="locus_id", how="left")

# Fill NA values if 'locus_id' from df1 is not found in df2
merged_df.fillna('NA', inplace=True)

# Save the merged data to a new file
merged_df.to_csv("8_pred_fisher_merged_dataset.txt", index=False, sep="\t")

print("Step 11 completed: predector results added'")


## Step 12: add known effectors--------------------------------------------------------------
#import csv
#import os
# Initialize new_data list with headers including the new 'known_effector' column
new_data = []

# Read known_effectors from file into a dictionary
known_effectors_dict = {}
if os.path.exists("input_files/known_effector.txt"): 
    with open("input_files/known_effector.txt", "r") as ke_file:
        for line in ke_file:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                known_effectors_dict[parts[0]] = parts[1]

# Process the dataset and add known_effector information
with open("8_pred_fisher_merged_dataset.txt", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    header = next(reader)
    new_data.append(header + ["known_effector"])  # Add the new column header
    
    for row in reader:
        # Retrieve the known_effector value based on locus_id, default to empty string if not found
        known_effector = known_effectors_dict.get(row[header.index("locus_id")], "")
        row.append(known_effector)
        new_data.append(row)

# Write the new data, including the known_effector column, back to a new file
with open("9_complete_isoform_list.txt", "w", newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerows(new_data)

print("Step 12 completed. The dataset with known effector added has been saved.")



## Step 13: final locus list--------------------------------------------------------------------------------------------
#import pandas as pd

# Load the data into a pandas DataFrame
data = pd.read_csv("9_complete_isoform_list.txt", sep="\t") #####change data set here############

# Sort the data based on PF_score
sorted_data = data.sort_values(by="p-value_lowest", ascending=True)
# Remove rows with duplicate values in 'variable' column, while keeping the first occurrence
data_no_duplicates = sorted_data.drop_duplicates(subset='locus_id', keep='first')

# Sort the data based on PF_score
sorted_data_for_ranking = data_no_duplicates.sort_values(by="effector_score", ascending=False)

sorted_data_for_ranking.to_csv('10_complete_locus_list.txt', sep='\t', index=False) #####change data set here############



## step 14 - filtering
# Initialize new_data to store rows that pass the filters
new_data = []

# Process file and filter data based on thresholds
with open("10_complete_locus_list.txt", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    header = next(reader)  # Read the header
    new_data.append(header)  # Add the header to new_data

    # Find the indices for the columns you'll be checking
    p_value_idx = header.index("p-value_lowest")
    pred_score_idx = header.index("effector_score")
    cyst_idx = header.index("aa_c_number")
    total_aa_idx = header.index("residue_number")

    # Iterate over each row in the reader
    for row in reader:
        try:
            # Convert string to float for comparison
            p_value = float(row[p_value_idx])
            pred_score = float(row[pred_score_idx])
            cyst = float(row[cyst_idx])
            total_aa = float(row[total_aa_idx])
            
            # Apply your filtering criteria
            if cyst >= args.cyst and total_aa <= args.total_aa and pred_score >= args.pred_score and p_value <= args.p_value:
                # If the row passes all filters, add it to new_data
                new_data.append(row)
        except ValueError:
            # Skip rows with conversion errors (e.g., empty strings or non-numeric values)
            continue
        
# Write the processed data to a new file
with open("11_filtered_candidate_list.txt", "w", newline='') as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerows(new_data)

print("Dataset filtered based on provided thresholds is saved.")



## step 15: ranking of known effector after filtering
import pandas as pd

# Load the data
data_path = "11_filtered_candidate_list.txt"
data = pd.read_csv(data_path, sep="\t")

# Ensure that 'known_effector' column is treated as a string
data['known_effector'] = data['known_effector'].astype(str)

# Sort the data by effector_score in descending order
sorted_data = data.sort_values(by="effector_score", ascending=False)

# Add a column for ranking based on the sorting
sorted_data['Ranking'] = range(1, len(sorted_data) + 1)

# Filter rows where "known_effector" is not empty (non-null and not an empty string)
# Also check that 'known_effector' is not 'nan' after conversion to string
known_effectors = sorted_data[
    (sorted_data["known_effector"].str.strip() != "") &
    (sorted_data["known_effector"].str.lower() != 'nan')
]

# If no known effectors are present, return or print a message
if known_effectors.empty:
    print("No known effectors found.")
else:
    # Select and reorder the columns in the specified sequence
    final_columns = ['Ranking', 'known_effector', 'locus_id', 'effector_score', 'p-value_lowest']
    known_effectors_final = known_effectors[final_columns]

    # Print or save to file
    print(known_effectors_final)

    # Optional: Save the processed data to a new CSV file
    output_path = "12_known_effectors_ranking.txt"
    known_effectors_final.to_csv(output_path, sep="\t", index=False)
