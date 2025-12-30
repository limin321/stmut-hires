import os
import pandas as pd
import numpy as np
import logging

from stmut_hires.call_cnv.cnvscore_permutation import CNVScore
from stmut_hires.call_cnv.cnvscore_fdr import CNVScoreFDR
from stmut_hires.call_cnv.save_cnv import CNVWriter

logger = logging.getLogger(__name__)

class CNVSortedByBulkCNV:
    """Sort Spatial CNV by corresponding bulkCNV info."""
    def __init__(self, bulk_csv: str, pmtimes: int, output_dir: str, clusterSortedcnv: pd.DataFrame):
        self.bulk_csv = bulk_csv
        self.pmtimes = pmtimes
        self.output_dir = output_dir
        self.clusterSortedcnv = clusterSortedcnv
        self.bulk = None
        self.sorted_cdt = None
        self.caseCNVscore_processed = None
        self.permuts_long = None
        self.caseCNVscore_FDR = None

    def read_bulk_cnv(self) -> None:
        if self.bulk_csv is not None:
            try:
                bulk = pd.read_csv(self.bulk_csv)
                bulk['gainloss'] = bulk['gainloss'].astype(int)
                self.bulk = bulk
                logger.info(f"Bulk CNV file loaded successfully from {self.bulk_csv}")
            except FileNotFoundError:
                logger.error(f"Bulk CNV file not found at {self.bulk_csv}")
        else:
            logger.warning("Warning: No bulk CNV path provided.")

    def case_cnvscore_permutation(self):
        """Calculate Case CNVScore, Run permutation, FDR"""
        gene_counts = CNVScore.count_arm_genes(self.clusterSortedcnv)
        sorted_cdt, caseCNVscore,permuts_df = CNVScore.sort_cnv_score(self.clusterSortedcnv, self.bulk, gene_counts, self.pmtimes)
        # Case CNVscore & permutation
        caseCNVscore_processed = CNVScoreFDR.calculate_observed_fractions(caseCNVscore)
        permuts_long = CNVScoreFDR.calculate_permuted_fractions(permuts_df)
        caseCNVscore_FDR = CNVScoreFDR.calculate_observed_FDR(permuts_long, caseCNVscore_processed)
        self.sorted_cdt = sorted_cdt
        self.caseCNVscore_processed = caseCNVscore_processed
        self.permuts_long = permuts_long
        self.caseCNVscore_FDR = caseCNVscore_FDR
        return permuts_long, caseCNVscore_FDR

    def run_pipeline(self):
        """Run bulk sorted CNV pipeline"""
        self.read_bulk_cnv()
        if self.bulk is not None:
            self.case_cnvscore_permutation()
            CNVWriter.save_bulk_sorted_data(
                self.sorted_cdt,
                self.caseCNVscore_processed,
                self.permuts_long,
                self.caseCNVscore_FDR,
                self.output_dir
            )
            logger.info("Bulk sorted CNV pipeline execution complete.")
            return self.permuts_long, self.caseCNVscore_FDR
        else:
            logger.warning("CNV sort by bulkCNV pipeline skipped due to missing bulk CSV path.")
            return None, None




