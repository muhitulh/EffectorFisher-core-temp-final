"""
processPredector - A utility class for processing Predector input

This module merges in-memory predector with the resulst from EffectorFisher
"""

import os
import pandas as pd
import logging

# Logger setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class ProcessPredector:

    def __init__(self, input_dir: str = '00_input_files'):
        self.input_dir = input_dir
        self.predector_df = pd.DataFrame()
        self.merged_df = pd.DataFrame()

    def load_data_predector(self):
        """Load predector data from default path."""
        predector_file = os.path.join(self.input_dir, "0_predector_results.txt")
        try:
            self.predector_df = pd.read_csv(predector_file, sep="\t")
            logger.info(f"Loaded Predector file: {predector_file}")
        except Exception as e:
            logger.error(f"Failed to load Predector file: {e}")
            raise

    def merge_with_fisher(self, fisher_path: str):
        """Merge predector data with Fisher results based on 'locus_id'."""
        try:
            self.fisher_df = pd.read_csv(fisher_path, sep="\t")
            logger.info(f"Loaded Fisher file: {fisher_path}")
        except Exception as e:
            logger.error(f"Failed to load Fisher file: {e}")
            raise

        # Ensure locus_id is the first column
        self.fisher_df = self._move_locus_id_first(self.fisher_df)
        self.predector_df = self._move_locus_id_first(self.predector_df)

        self.merged_df = self.fisher_df.merge(self.predector_df, on="locus_id", how="left")
        self.merged_df.fillna('NA', inplace=True)
        logger.info("Step 10 completed: Predector results merged with Fisher output.")

    def _move_locus_id_first(self, df: pd.DataFrame) -> pd.DataFrame:
        """Move 'locus_id' column to the first position."""
        if 'locus_id' in df.columns:
            cols = ['locus_id'] + [col for col in df.columns if col != 'locus_id']
            return df[cols]
        return df

    def save_processed_data(self, output_path: str):
        """Save the merged dataframe to a specified output file."""
        if self.merged_df.empty:
            logger.warning("No data to save. Merge the data first.")
            return
        try:
            self.merged_df.to_csv(output_path, sep="\t", index=False)
            logger.info(f"Merged dataset saved to: {output_path}")
        except Exception as e:
            logger.error(f"Failed to save output: {e}")
            raise
