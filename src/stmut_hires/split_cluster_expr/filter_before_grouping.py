import pandas as pd
import numpy as np
from scipy.sparse import csc_matrix
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import logging

class FilterOutCrappyData:
    """ 
    Filter out crappy data before grouping to avoid fan-shape issue.
    Low quality data grouping together will lead to fan-shape when visualizing spatially.
    mat_dict is the dictionary from h5_reader.py
    Either automatically generate filter-cutoff or user provide the cutoff
    """
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)

    def get_filter_cutoff(self, mat_dict, bw_method=0.1, filter_cutoff=None):
        """
        Using Kernel Density Estimation (KDE) Valley Detection to automatically select a cutoff to filter-out gene
        expression h5 matrix.
        bw_method: default value is 0.1, increasing bandwidth value, more data will be filter-out. 
        
        """
        barcodes = mat_dict['barcodes']
        indptr = mat_dict['indptr']
        nnz_per_barcode = np.diff(indptr)
        df = pd.DataFrame({
            "barcodes":barcodes,
            "gene_counts":nnz_per_barcode
        })

        counts = df['gene_counts']

        if filter_cutoff is None:
            self.logger.info("Use kde predict cutoff")
            # 1. Increase bandwidth (e.g., 0.3 or 0.4) to smooth over micro-fluctuations near 0
            kde = gaussian_kde(counts, bw_method=bw_method)
            x_axis = np.linspace(0, max(counts), 1000)
            kde_values = kde(x_axis)
            
            valleys, _ = find_peaks(-kde_values) # Find local minima (valleys)

            if len(valleys) > 0: # Safely pick the valley that is furthest along or the first true macro-valley
                # Pick the first valley detected after the smooth curve drops
                cutoff = x_axis[valleys[0]]
            else: # Fallback if no valley is found after heavy smoothing
                cutoff = np.percentile(counts, 5)

        else: # user-provided filter cutoff
            self.logger.info("Use user provided cutoff")
            cutoff = filter_cutoff
        
        self.logger.info(f"Barcodes have less than {cutoff} genes are filtered out.")

        return df, cutoff

    def filtered_dict(self, mat_dict, cutoff, df):
        "Filter filtered_feature_bc_matrix.h5 based on user-provided cutoff"
        mat = csc_matrix(
            (
                mat_dict['data'],
                mat_dict['indices'],
                mat_dict['indptr']
            ),
            shape=tuple(mat_dict['shape'])
        )

        nnz = np.diff(mat_dict['indptr'])
        keep_mask = nnz >= cutoff
        filtered_mat = mat[:, keep_mask]
        filtered_barcodes = mat_dict['barcodes'][keep_mask]
        filtered_dict = {
            'data': filtered_mat.data,
            'indices': filtered_mat.indices,
            'indptr': filtered_mat.indptr,
            'shape': filtered_mat.shape,
            'barcodes': filtered_barcodes,
            'genes': mat_dict['genes']
        }
        
        filtered_df = df[df["gene_counts"] >= cutoff]
        return filtered_dict, filtered_df

