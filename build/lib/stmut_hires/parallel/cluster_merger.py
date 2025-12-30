import pyarrow.parquet as pq
import logging

from stmut_hires.cluster_processor.merge_cluster_barcodes import BarcodeMerger

class ClusterMerger:
    """Responsible for processing individual clusters"""
    def __init__(self, coords, txt_output_dir, ensembl_path, cutoff=1000, window=100):
        self.coords = coords
        self.txt_output_dir = txt_output_dir
        self.ensembl_path = ensembl_path
        self.cutoff = cutoff
        self.window = window
        self.logger = logging.getLogger(__name__)

    def process_cluster(self, cluster_path):
        """Merge barcodes for each cluster"""
        self.logger.info(f"Working on {cluster_path}")
        # Extract all barcodes of the cluster 
        table = pq.read_table(cluster_path)  # loads no rows
        barcodes = table.schema.names

        # Subset the spatial barcodes of the cluster
        cluster_coor = self.coords[self.coords['barcode'].isin(barcodes)]
                
        merger = BarcodeMerger(
            cluster_spatial_coords=cluster_coor, 
            expression_file=cluster_path, 
            ensembl_file=self.ensembl_path, 
            output_dir=self.txt_output_dir, 
            cutoff=self.cutoff, 
            window=self.window
        )
        merger.merge_barcodes()

