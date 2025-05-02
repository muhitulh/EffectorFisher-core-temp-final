#!/usr/bin/env python3
"""
EffectorFisher Core Script

This script processes phenotype data for the EffectorFisher pipeline.
It can handle either quantitative or qualitative data.
"""

import os
import argparse
import pandas as pd
from utils.phenotypeProcessor import PhenotypeProcessor


def parse_arguments():
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(description='Process phenotype data for EffectorFisher')
    
    parser.add_argument(
        '--data-type',
        choices=['quantitative', 'qualitative'],
        default='quantitative',
        help='Type of data to process (quantitative or qualitative)'
    )
    
    parser.add_argument(
        '--input-dir',
        default='00_input_files',
        help='Directory containing input files'
    )
    
    parser.add_argument(
        '--output-dir',
        default='output',
        help='Directory for output files'
    )
    
    parser.add_argument(
        '--save',
        action='store_true',
        help='Save processed data to files'
    )
    
    return parser.parse_args()


def main():
    """
    Main function to run the EffectorFisher core process.
    """
    # Parse command line arguments
    args = parse_arguments()
    
    print(f"Processing {args.data_type} phenotype data...")
    
    # Initialize the phenotype processor
    processor = PhenotypeProcessor(input_dir=args.input_dir)
    
    try:
        # Step 1: Process the data (load, convert, split)
        trait_dataframes = processor.process_data(data_type=args.data_type)
        
        # Step 2: Report processed trait summary
        print(f"Successfully processed {len(trait_dataframes)} traits:")
        for trait_name, df in trait_dataframes.items():
            print(f"  - {trait_name}: {len(df)} samples (before cleaning)")

        # Step 3: Optionally clean & save
        if args.save:
            print(f"Cleaning and saving processed data to '{args.output_dir}'...")
            cleaned_outputs = processor.save_processed_data(output_dir=args.output_dir)

            for fname, df in cleaned_outputs.items():
                print(f"  - {fname}: {len(df)} samples (after cleaning)")

            print("Data cleaned and saved successfully.")
    
    except Exception as e:
        print(f"Error: {str(e)}")
        return 1
    
    print("Processing complete.")
    return 0


if __name__ == "__main__":
    exit(main())

