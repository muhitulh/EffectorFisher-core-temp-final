#!/usr/bin/env python3
"""
EffectorFisher Core Script

This script processes phenotype and variant data for the EffectorFisher pipeline.
It supports both quantitative and qualitative phenotype data and variant filtering.
"""

import os
import argparse
import pandas as pd
from utils.phenotypeProcessor import PhenotypeProcessor
from utils.variantProcessor import VariantProcessor


def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(description='Process phenotype and variant data for EffectorFisher')

    parser.add_argument(
        '--data-type',
        choices=['quantitative', 'qualitative'],
        default='quantitative',
        help='Type of phenotype data to process'
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
        '--min-variant',
        type=int,
        default=5,
        help='Minimum variant frequency for filtering'
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
    args = parse_arguments()

    # ─────── PHENOTYPE PROCESSING ───────
    print(f"Processing {args.data_type} phenotype data...")
    phenotype_processor = PhenotypeProcessor(input_dir=args.input_dir)

    try:
        trait_dataframes = phenotype_processor.process_data(data_type=args.data_type)
        print(f"Successfully processed {len(trait_dataframes)} traits:")
        for trait_name, df in trait_dataframes.items():
            print(f"  - {trait_name}: {len(df)} samples (before cleaning)")
    except Exception as e:
        print(f"[Phenotype Error] {str(e)}")
        return 1

    # ─────── VARIANT PROCESSING ───────
    print(f"\nProcessing variant data with min_variant = {args.min_variant}...")
    variant_processed = VariantProcessor(input_dir=args.input_dir)

    try:
        variant_processed.load_data()
        variant_processed.filter_by_variant_frequency(min_var=args.min_variant)
        print(variant_processed.filtered_df)
    except Exception as e:
        print(f"[Variant Error] {str(e)}")
        return 1

    # ─────── SAVE ALL OUTPUTS ───────
    if args.save:
        print(f"\nSaving all processed data to '{args.output_dir}'...")

        # Save phenotype
        try:
            phenotype_outputs = phenotype_processor.save_processed_data(output_dir=args.output_dir)
            for fname, df in phenotype_outputs.items():
                print(f"  - {fname}: {len(df)} samples (phenotype)")
        except Exception as e:
            print(f"[Save Error - Phenotype] {str(e)}")
            return 1

        # Save variant
        try:
            variant_output_file = os.path.join(args.output_dir, '0_filtered_combined_variants.txt')
            variant_processor.save_processed_data(output_file=variant_output_file)
            print(f"  - {variant_output_file}: {len(variant_processor.filtered_df)} samples (variants)")
        except Exception as e:
            print(f"[Save Error - Variant] {str(e)}")
            return 1

        print("All data cleaned and saved successfully.")

    print("\nProcessing complete.")
    return 0


if __name__ == "__main__":
    exit(main())
