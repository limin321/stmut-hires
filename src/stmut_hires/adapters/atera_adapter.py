import zarr
import pandas as pd
import numpy as np
import os

class AteraAdapter:
    """ 
    Convert Atera data format to be compatible with the main pipeline - Visiumhd format
    Requires 3 files:
    cell_feature_matrix.h5 - no need to change. same structure as filtered_feature_barcodes.h5
    cells.parquet, analysis.zarr.zip -- these two are required to get cluster_file (Graph-based.csv), and spatial_file (tissue_position_lsit.csv)
    """
    def __init__(self, output_dir, dry_run=False):
        self.output_dir = output_dir
        self.dry_run = dry_run
        if not dry_run:
            full_path = os.path.join(output_dir, "atera")
            os.makedirs(full_path, exist_ok=True)

    """Convert Atera output to VisiumHD Compatible Format """
    def get_graph_based(self, analysis_zarr, cells_parquet):
        cluster_file = f"{self.output_dir}/atera/Graph-based.csv"
        if self.dry_run:
            return cluster_file

        cells_df = pd.read_parquet(cells_parquet)
        barcodes = cells_df['cell_id'].values
        store = zarr.storage.ZipStore(analysis_zarr, mode='r')
        root = zarr.open(store)
        indices = root['cell_groups/0/indices'][:]
        indptr = root['cell_groups/0/indptr'][:]
        cluster_assignments = np.zeros(len(barcodes), dtype=int)

        for cluster_idx in range(len(indptr) - 1):
            start = indptr[cluster_idx]
            end = indptr[cluster_idx + 1]
            affected_cell_indices = indices[start:end]
            # Assign cluster ID (using 1-based indexing to match Explorer usually)
            cluster_assignments[affected_cell_indices] = cluster_idx + 1

        result_df = pd.DataFrame({
            'Barcode': barcodes,
            'Graph-based': cluster_assignments
        })

        result_df['Graph-based'] = "Cluster " + result_df['Graph-based'].astype(str)
        result_df.to_csv(cluster_file, index=False)
        return cluster_file


    def get_tissue_position_parquet(self, cells_parquet):
        spatial_file = f"{self.output_dir}/atera/tissue_positions.parquet"
        if self.dry_run:
            return spatial_file
            
        atera = pd.read_parquet(cells_parquet)
        atera_sub = atera[["cell_id", "x_centroid", "y_centroid"]]
        atera_sub = atera_sub.rename(
            columns={
                'cell_id': "barcode",
                'x_centroid': "array_col",
                'y_centroid': "array_row"
            }
        )
        atera_sub.to_parquet(spatial_file)
        return spatial_file