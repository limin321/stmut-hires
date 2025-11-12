import numpy as np
from sklearn.neighbors import NearestNeighbors
from tqdm.auto import tqdm # For automatic selection of notebook/console bar

class BarcodeSorter:
    "Responsible for updating window-size for grouping"

    @staticmethod
    def dynamic_update_window(barcodes, window, cutoff, coords, current_idx, parquet_file):
        # Initialize nbrs outside the loop
        nbrs = NearestNeighbors(n_neighbors=min(window, len(coords)), algorithm="kd_tree").fit(coords)
        
        # Initial calculation before the loop
        distances, indices = nbrs.kneighbors([coords[current_idx]])
        sorted_indices = indices[0][np.argsort(distances[0])]
        sorted_barcodes = barcodes[sorted_indices]

        # subset expression data and calculate non-zero gene count
        sub_expr = parquet_file.read(columns=sorted_barcodes).to_pandas()
        sub_expr_sum = sub_expr.sum(axis=1)
        sub_total_non_zero_count = (sub_expr_sum != 0).sum()

        # Loop to dynamically update the window and re-evaluate
        while sub_total_non_zero_count < cutoff and window <= len(coords):
            tqdm.write("Dynamically update window size for this cluster.")
            # Dynamic window update logic
            window *= 2  # Double the window size
            tqdm.write(f"New window size is {window}")
            # Re-fit the model with the larger window for the next iteration
            nbrs = NearestNeighbors(n_neighbors=min(window, len(coords)), algorithm="kd_tree").fit(coords)
            
            # Get nearest neighbors with the new window size
            distances, indices = nbrs.kneighbors([coords[current_idx]])
            sorted_indices = indices[0][np.argsort(distances[0])]
            sorted_barcodes = barcodes[sorted_indices]

            # subset expression data and calculate non-zero gene count
            sub_expr = parquet_file.read(columns=sorted_barcodes).to_pandas()
            sub_expr_sum = sub_expr.sum(axis=1)
            sub_total_non_zero_count = (sub_expr_sum != 0).sum()

        return sorted_barcodes

