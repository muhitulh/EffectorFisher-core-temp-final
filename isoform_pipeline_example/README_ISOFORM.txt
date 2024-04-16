Input file in directory - 01_input_from_metaeuk:
- Gene-specific protein FASTA file from metaeuk screening using pan-gene set
- Example: 01_input_from_metaeuk/Effector_screen_SNOO_001050A_SN15_SNOG_001050A_prot.fas

Pipeline scripts:

Script 1: 1_isoform_run.sh
   - This script runs the Perl script "cluster_exact_seqs_fasta.pl" on multiple input files to extract isoforms.
   - It generates multiple files, including *_prot.fas.byseqid_table, which is an isoform frequency table.

Script 2: 2_rename_extract_seq_and_combine_isoforms.py
	section 1:
   - This script changes the isoform sequence headers to isoform IDs.
   - It generates a new CSV file with renamed isoform IDs and their corresponding sequences.
   - This allows for easy lookup of the sequence for a particular isoform later on.
	
	section 2:
   - It combines all the gene-specific isoform tables into a combined master isoform table.
   - The output file is located at 02_isoform_combined_table/combined_isoform_output.txt.
   - This combined table will be used as input for effectorfisher.

Script 3: 3_complete_isoform_id_with_corresponding_seq.sh
   - This script generates a new CSV file with renamed isoform IDs and their corresponding sequences.
   - The output file is located at 03_isoform_seq/isoform_seq.txt.
   - This file allows for easy lookup of the sequence for a particular isoform later on.
