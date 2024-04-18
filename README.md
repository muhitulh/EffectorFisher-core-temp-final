# EffectorFisher_test

Usage:Run the script with the following command:
python effectorfisher_pipleline.py --data_type quantitative --min_iso 5 --cyst 2 --pred_score 2 --total_aa 300 --p_value 0.05

Options:
Must to include:<br>
--data_type qualitative/quantitative [you have to select what kind of phenotypic data you have, see example in input_files]

Important:<br>
--min_iso [minimum isoform number] (default = 5)

Optional:<br>
--cyst [Cysteine count threshold] (default = 2)
--pred_score [Prediction score threshold] (default = 2)
--total_aa [Total amino acid count threshold] (default = 300)
--p_value [P-value threshold] (default = 0.05)
If you don't define these thresholds, the script will use the default values.

##Output<br>
#Main output<br>
9_complete_isoform_list.txt
10_complete_locus_list.txt

#Additional output<br>
11_filtered_candidate_list.txt: Based on the default or specified filters. Alternatively, you can apply filters to 10_complete_locus_list.txt as required.
12_known_effectors_ranking.txt: Contains the ranking of known effectors if you provide a known effector input file.

##Input files<br>
0_combined_isoform.txt: This file will be found after running the isoform pipeline in the directory EffectorFisher_test/00_Isoform_table_pipeline_example/02_isoform_combined_table/0_combined_isoform.txt.

0_phenotype_data.txt: You need to prepare this file as shown in the example.

0_predector_results.txt: This file will be found after running the predector pipeline.
known_effector.txt: You can provide known effector IDs and names as shown in the example. If not provided, it will be empty.
