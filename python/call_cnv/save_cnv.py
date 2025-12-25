import pandas as pd
import os 
import logging
from typing import Optional

logger = logging.getLogger(__name__)

class CNVWriter:
    """Responsible for saving CNV-related files."""
    @staticmethod
    def _ensure_tables_dir(output_dir: str) -> str:
        """Ensure the `tables` directory exists under the output directory."""
        tables_dir = os.path.join(output_dir, "tables")
        os.makedirs(tables_dir, exist_ok=True)
        return tables_dir

    @staticmethod
    def save_merged_df(merged_df: Optional[pd.DataFrame], output_dir: str) -> None:
        """
        Save the merged barcodes summary DataFrame to a CSV file.

        Args:
            merged_df (Optional[pd.DataFrame]): DataFrame containing merged barcode summaries.
            output_dir (str): Directory where the CSV file will be saved.

        Notes:
            - Saves the CSV under <output_dir>/tables/cluster_barcodes_summary.csv.
            - Logs success or warning if DataFrame is None.
        """
        if merged_df is None:
            logger.warning("No merged_df to save.")
            
        tables_dir = CNVWriter._ensure_tables_dir(output_dir)
        save_path = os.path.join(tables_dir, "cluster_barcodes_summary.csv")
        try:
            merged_df.to_csv(save_path, index=False)
            logger.info(f"Saved merged summary to {save_path}")
        except Exception as e:
            logger.error(f"Failed to save merged_df: {e}")
            raise

    @staticmethod
    def save_cluster_sorted_cnv(cdt_sorted, cdt_meta, output_dir):
        """
        Save CNV data organized by gene expression clusters and total reads.

        Args:
            cdt_sorted (pd.DataFrame): Cluster-sorted CNV data including 'cluster', 'annotate', 'TotalRDs', 'barcode'.
            cdt_meta (pd.DataFrame): Metadata corresponding to CNV matrix.
            output_dir (str): Directory to save output files.

        Returns:
            pd.DataFrame: Merged dataframe combining metadata and CNV data.
        """
        tables_dir = CNVWriter._ensure_tables_dir(output_dir)
        output_file = os.path.join(tables_dir, "CNVs_OrganizedByGEcluster_UMIcount.cdt")
        try:
            cdt_sorted = cdt_sorted.drop(['cluster', 'annotate', 'TotalRDs'], axis=1)
            cdt_sorted_indexed = cdt_sorted.set_index('barcode').transpose()
            merged = pd.concat([cdt_meta.reset_index(drop=True), cdt_sorted_indexed.reset_index(drop=True)], axis=1)
            merged.to_csv(output_file, sep = "\t", index = False)
            logger.info(f"Cluster-sorted CNV saved to {output_file}")
            return merged
        except Exception as e:
            logger.error(f"Failed to save cluster-sorted CNV: {e}")
            raise

    @staticmethod
    def save_qt(all_barcodes_quintiles,output_dir):
        """
        Save barcode quintile data.

        Args:
            all_barcodes_quintiles (pd.DataFrame): DataFrame containing CNV quintile assignments.
            output_dir (str): Directory to save output files.
        """
        tables_dir = CNVWriter._ensure_tables_dir(output_dir)
        output_file = os.path.join(tables_dir, "CNVs_RankedbySimilaritytoDNA_Quintiles4Loupe.csv")
        try:
            output_dir_tables = os.path.join(output_dir, "tables")
            all_barcodes_quintiles.to_csv(output_file, index=False)
            logger.info(f"Saved barcode quintiles to {output_file}")
        except Exception as e:
            logger.error(f"Failed to save quintile data: {e}")
            raise
    
    @staticmethod
    def save_bulk_sorted_data(
        sorted_cdt: Optional[pd.DataFrame],
        caseCNVscore_processed: Optional[pd.DataFrame],
        permuts_long: Optional[pd.DataFrame],
        caseCNVscore_FDR: Optional[pd.DataFrame],
        output_dir: str
    ) -> None:
        """
        Save bulk CNV and permutation data.

        Args:
            sorted_cdt (pd.DataFrame): Bulk sorted CNV matrix.
            caseCNVscore_processed (pd.DataFrame): Processed CNV scores for histogram plotting.
            permuts_long (pd.DataFrame): Permutation results.
            caseCNVscore_FDR (pd.DataFrame): FDR-adjusted CNV scores.
            output_dir (str): Directory to save output files.
        """
        if sorted_cdt is None:
            logger.warning("Bulk sorted data not generated yet. Run case_cnvscore_permutation() first.")
            return
        tables_dir = CNVWriter._ensure_tables_dir(output_dir)
        try:
            sorted_cdt.to_parquet(os.path.join(tables_dir, "CNVs_RankedBySimilarityToDNA.cdt.parquet"))
            caseCNVscore_processed.to_csv(os.path.join(tables_dir, "CNVs_RankedBySimilarityToDNA_CNVscoreHistogram.csv"), index=False)
            permuts_long.to_parquet(os.path.join(tables_dir, "permut_CNVscores.parquet"))
            caseCNVscore_FDR.to_parquet(os.path.join(tables_dir, "caseCNVScore.parquet"))
            logger.info("Bulk sorted CNV data and permutations saved successfully.")
        except Exception as e:
            logger.error(f"Failed to save bulk sorted CNV data: {e}")
            raise



