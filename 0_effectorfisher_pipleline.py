## for thresholds set up in command line
#usage: python 0_effectorfisher_pipleline.py --min_iso 5 --cyst 2 --pred_score 2 --total_aa 300 --p_value 0.05


import pandas as pd
import argparse
import glob
import os

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




## Step 1 &2: create cultivar specific files--------------------------original step 2 deleted, as that included here

#import pandas as pd

# Load the file, ensuring that the first column is treated as the 'ID' column
# If the file does not have headers, uncomment the following line and comment out the pd.read_csv line below it
# df = pd.read_csv('0_phenotype_data.txt', sep='\t', header=None, names=['ID', 'Cultivar1', 'Cultivar2', ...])
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

print("Cultivar-specific files with 'disease' column have been created.")





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
filtered_df.to_csv('0_filtered_combined_isoform.txt', index=False, sep='\t')





####ste 5: ## step 5: merge cultivar-specific file with combine isoform file by ID

#import pandas as pd
#import glob
#import os


# File to be concatenated
additional_file = '0_filtered_combined_isoform.txt'

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

    print(f"Concatenated file saved: {concatenated_file_path}")
    
    
    
    

## Step 6: change 0 = A and 1 = P in those files----------------------------------------------
#import pandas as pd
#import glob

# Using glob to find all files starting with '1_data' in the specified directory
data_files = glob.glob('2_merged_data*.txt')

# Processing each file
for data_file in data_files:
    data_df = pd.read_csv(data_file, sep='\t')

    # Looping through all columns after the second one
    for col in data_df.columns[2:]:
        # Replace 1 with 'P' and 0 with 'A'
        data_df[col] = data_df[col].apply(lambda x: 'P' if x == 1.0 else ('A' if x == 0.0 else x))

    # Save the modified dataframe back to the file
    data_df.to_csv(data_file, sep='\t', index=False)
    
    
    
    
    
    
## step 7: contingency table------------------------------------------
#import pandas as pd

for i in range(1, 3):  # Loop from 1 to 12 for your 12 data files----------------------------------neeed to be automated here 
    # Dynamically create file name
    file_name = f'2_merged_data{i}.txt'

    # Load data
    data = pd.read_csv(file_name, delimiter='\t')  # Adjust delimiter if needed

    # Melt the data
    melted = pd.melt(data, id_vars=['disease'], var_name='isoform', value_name='value')

    # Pivot the table to get the contingency table
    pivot_data = melted.pivot_table(index='isoform', columns=['disease', 'value'], aggfunc='size', fill_value=0)

    # Create new headers
    pivot_data.columns = ['{}-{}'.format(disease, value) for disease, value in pivot_data.columns]

    # Filter columns
    desired_columns = ['high-A', 'high-P', 'low-A', 'low-P']
    pivot_data = pivot_data[desired_columns]

    # Remove the first row (which contains 'ID')
    pivot_data = pivot_data.drop(pivot_data.index[0])

    # Write the contingency table to a file
    contingency_table_file = f'3_contingency_table_data{i}.txt'
    pivot_data.to_csv(contingency_table_file, sep='\t', index=True)

    # Rename the columns
    column_rename = {'high-A': 'c', 'high-P': 'a', 'low-A': 'd', 'low-P': 'b'}
    pivot_data.rename(columns=column_rename, inplace=True)

    # Write the reformatted contingency table to a new file
    hypergeo_data_file = f'4_hypergeo_data{i}.txt'
    pivot_data.to_csv(hypergeo_data_file, sep='\t', index=True)

    print(f"Script completed for dataset {i}. The formatted contingency table has been written to '{contingency_table_file}' and '{hypergeo_data_file}'")
    
    
    
    
## Step 8: hypergeomatric test-------------------------------------------------
    
### automated: Hypergeo test - loop over many datasets together [automated]
###loop over all dataset
#import math

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

# Determine the range of your datasets
start_dataset_number = 1  # start number -------------------------------------------------------need to change here (need to be automated)
end_dataset_number = 2   # end number ---------------------------------------------------------need to change here (need to be automated)

# Pre-calculate factorials
fact = calculate_factorials(1000000)

# Loop over the datasets
for i in range(start_dataset_number, end_dataset_number + 1):
    input_file = f'4_hypergeo_data{i}.txt'
    output_file = f'5_output_hypergeo_data{i}.txt'
    print(f"Processing dataset number {i}...")
    process_dataset(input_file, output_file, fact)

print("All datasets processed.")





#step 9 & 10 merge all the cultivar-specific p-value (master data) and add lowest p-value col-----------------------------------------------------

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

print("Script completed. Merged data with lowest p-value written to '7_merged_lowest_p-value.txt'")




##step 11: ## Step 11: need to create col locus_id from isoform_id---------------------------------------
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





###step 12: combine fisher result with predector result----------------------------------------(try with raw predector output)
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

print("Merging complete!")





## Step 13: calculate predfish value--------------------------------------------------------------

#import csv
#import argparse
#import os

# Setup argument parser
## you can change default values here---------------------------------------------------------------------------------------------
#parser = argparse.ArgumentParser(description="Filter dataset based on provided thresholds.")
#parser.add_argument('--cyst', type=float, default=2, help='Cysteine count threshold')
#parser.add_argument('--pred-score', type=float, default=2, help='Prediction score threshold')
#parser.add_argument('--total-aa', type=float, default=300, help='Total amino acid count threshold')
#parser.add_argument('--p-value', type=float, default=0.05, help='P-value threshold')

# Parse arguments
#args = parser.parse_args()

# Initialize variable to find max pred_score
max_pred_score = 0

# Read known_effectors from file into a dictionary
known_effectors_dict = {}
if os.path.exists("input_files/known_effector.txt"): ##---------------------------------------------if you have different known_effector file, need to change there 
    with open("input_files/known_effector.txt", "r") as ke_file:
        for line in ke_file:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                known_effectors_dict[parts[0]] = parts[1]


# Read the file once to find the maximum prediction score
with open("8_pred_fisher_merged_dataset.txt", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    header = next(reader)
    pred_score_idx = header.index("effector_score")

    # Indexes based on header names
    variable_idx = header.index("locus_id")
    p_value_idx = header.index("p-value_lowest")
    pred_score_idx = header.index("effector_score")
    cyst_idx = header.index("aa_c_number")
    total_aa_idx = header.index("residue_number")

    ## To get the max predscore
    for row in reader:
        try:
            score = float(row[pred_score_idx])
            if score > max_pred_score:
                max_pred_score = score
        except ValueError:
            continue  # Skip rows with invalid data

# Re-open the file for processing rows based on thresholds
new_data = [header + ["convt_p-value", "nor_pred_score", "PF_score", "known_effector"]]

##unfiltered isoform output with all genes----------------------------------

with open("8_pred_fisher_merged_dataset.txt", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    next(reader)  # Skip header row

    for row in reader:
        try:
            # Calculate convt_p_value, nor_pred_score, and PF_score directly without the filtering conditions
            convt_p_value = 1 - float(row[p_value_idx])
            nor_pred_score = (float(row[pred_score_idx]) - args.pred_score) / (max_pred_score - args.pred_score)
            PF_score = (nor_pred_score * 0.5) + (convt_p_value * 0.5)

            # If known_effectors_dict is not empty, use it to set known_effector
            if known_effectors_dict:
                known_effector = known_effectors_dict.get(row[variable_idx], "")
            else:
                known_effector = ""

            row.extend([str(convt_p_value), str(nor_pred_score), str(PF_score), known_effector])
            new_data.append(row)
        except ValueError:
            continue  # Skip row if conversion fails

# Write the processed data to a new file
with open("9_pred_fisher_complete_isoform_set.txt", "w", newline='') as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerows(new_data)


print("Unfiltered dataset saved with all the parameters.")



#---------------------------

# Setup argument parser
## you can change default values here---------------------------------------------------------------------------------------------
#parser = argparse.ArgumentParser(description="Filter dataset based on provided thresholds.")
#parser.add_argument('--cyst', type=float, default=2, help='Cysteine count threshold')
#parser.add_argument('--pred_score', type=float, default=2, help='Prediction score threshold')
#parser.add_argument('--total_aa', type=float, default=300, help='Total amino acid count threshold')
#parser.add_argument('--p_value', type=float, default=0.05, help='P-value threshold')

# Parse arguments
#args = parser.parse_args()

# Initialize variable to find max pred_score
max_pred_score = 0

# Read known_effectors from file into a dictionary
known_effectors_dict = {}
if os.path.exists("known_effector.txt"): ##---------------------------------------------if you have different known_effector file, need to change there 
    with open("known_effector.txt", "r") as ke_file:
        for line in ke_file:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                known_effectors_dict[parts[0]] = parts[1]


# Read the file once to find the maximum prediction score
with open("8_pred_fisher_merged_dataset.txt", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    header = next(reader)
    pred_score_idx = header.index("effector_score")

    # Indexes based on header names
    variable_idx = header.index("locus_id")
    p_value_idx = header.index("p-value_lowest")
    pred_score_idx = header.index("effector_score")
    cyst_idx = header.index("aa_c_number")
    total_aa_idx = header.index("residue_number")

    ## To get the max predscore
    for row in reader:
        try:
            score = float(row[pred_score_idx])
            if score > max_pred_score:
                max_pred_score = score
        except ValueError:
            continue  # Skip rows with invalid data

# Re-open the file for processing rows based on thresholds
new_data = [header + ["convt_p-value", "nor_pred_score", "PF_score", "known_effector"]]

###filtered isoform output-----------------------------------------------------------------------------------------------
with open("8_pred_fisher_merged_dataset.txt", "r") as file:
    reader = csv.reader(file, delimiter="\t")
    next(reader)  # Skip header row

    for row in reader:
        try:
            if float(row[cyst_idx]) < args.cyst or float(row[total_aa_idx]) > args.total_aa or float(row[pred_score_idx]) < args.pred_score or float(row[p_value_idx]) > args.p_value:
                continue
            convt_p_value = 1 - float(row[p_value_idx])
            nor_pred_score = (float(row[pred_score_idx]) - args.pred_score) / (max_pred_score - args.pred_score)
            PF_score = (nor_pred_score * 0.5) + (convt_p_value * 0.5)

            # Look up known_effector using the locus_id from the row
            if known_effectors_dict:  # Check if known_effectors_dict is not empty
                known_effector = known_effectors_dict.get(row[variable_idx], "")
            else:
                known_effector = ""

            row.extend([str(convt_p_value), str(nor_pred_score), str(PF_score), known_effector])
            new_data.append(row)
        except ValueError:
            continue  # Skip row if conversion fails

# Write the processed data to a new file
with open("10_pred_fisher_filtered_isoform_set.txt", "w", newline='') as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerows(new_data)

print("dataset filtered based on provided thresholds is saved.")




## Step 14: final locus list--------------------------------------------------------------------------------------------
#import pandas as pd

# Load the data into a pandas DataFrame
data = pd.read_csv("9_pred_fisher_complete_isoform_set.txt", sep="\t") #####change data set here############
# Remove rows with duplicate values in 'variable' column, while keeping the first occurrence
data_no_duplicates = data.drop_duplicates(subset='locus_id', keep='first')
# Save the modified DataFrames back to files
data_no_duplicates.to_csv('9_pred_fisher_complete_gene_set.txt', sep='\t', index=False) #####change data set here############


# Load the data into a pandas DataFrame
data = pd.read_csv("10_pred_fisher_filtered_isoform_set.txt", sep="\t") #####change data set here############
# Remove rows with duplicate values in 'variable' column, while keeping the first occurrence
data_no_duplicates = data.drop_duplicates(subset='locus_id', keep='first')
# Save the modified DataFrames back to files
data_no_duplicates.to_csv('10_pred_fisher_filtered_gene_set.txt', sep='\t', index=False) #####change data set here############

print("Complete candidates and top ranked EffectorFisher candidates lists are saved")




## step 15: ranking by PF_score------------------------------------------------
###ranking based on three methods

#import pandas as pd

# Load the data into a pandas DataFrame
data = pd.read_csv("10_pred_fisher_filtered_gene_set.txt", sep="\t") #####change data set here############

# Sort the data based on PF_score
sorted_data = data.sort_values(by="PF_score", ascending=False)

# Find and populate the rows where "known_effector" column is not empty
effectors_table = [["Ranking", "Known_effector", "p-value_lowest", "effector_score", "PF_score"]]
for i, (_, row) in enumerate(sorted_data.iterrows(), start=1):  # Using iterrows to iterate over DataFrame rows
    if pd.notnull(row["known_effector"]):  # Using pandas' notnull to check for non-empty values
        effectors_table.append([i, row["known_effector"], row["p-value_lowest"], row["effector_score"],row["PF_score"]])

# Convert the effectors_table list to a DataFrame
effectors_df = pd.DataFrame(effectors_table[1:], columns=effectors_table[0])

# Save the modified DataFrames back to files
effectors_df.to_csv('11_ranking_known_effector_by_PF.txt', sep='\t', index=False) #####change data set here############

print("Final EffectorFisher ranking of known effector saved")

