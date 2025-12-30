from .cluster_reader import ClusterReader
from .h5_reader import H5MatrixReader
from .matrix_builder import ExpressionMatrixBuilder
from .gene_writer import GeneDataWriter
from tqdm import tqdm
import logging

""" 6. Main Processor(Orchestrator) """
class ClusterExpressionProcessor:
    """Orchestrate the cluster expression processing workflow """
    def __init__(self, config, dry_run=False):
        self.config = config
        self.dry_run = dry_run
        self.logger = logging.getLogger(__name__)
        self.cluster_reader = ClusterReader()
        self.h5_reader = H5MatrixReader(config.exp_h5)
        self.matrix_builder = ExpressionMatrixBuilder()
        self.gene_writer = GeneDataWriter(config.output_dir, dry_run=dry_run)

    def process(self):
        """ Main processing workflow """
        self.logger.info("Starting cluster expression processing")

        # Read cluster data
        self.logger.info("Reading cluster data...")
        cluster_data = self.cluster_reader.read_clusters(self.config.clusterf)
        self.logger.info(f"Found {len(cluster_data)} clusters")

        # Read hdf5 matrix
        self.logger.info("Reading HDF5 matrix data...")
        matrix_data = self.h5_reader.read_matrix_data()
        self.logger.info(f"Matrix shape: {matrix_data['shape']}")

        # Build expression Matrix
        self.logger.info("Building expression matrix...")
        expr_matrix = self.matrix_builder.build_sparse_matrix(matrix_data)
        barcode_index_map = self.matrix_builder.build_barcodes_index_map(matrix_data)
        self.logger.info(f"Built barcode index map with {len(barcode_index_map)} barcodes")

        # Save gene ids
        self.gene_writer.save_gene_ids(matrix_data['genes'])

        # Process each cluster
        self._process_clusters(
            cluster_data,
            barcode_index_map, 
            expr_matrix, 
            matrix_data['barcodes']
        )
        
        self.logger.info("Cluster expression processing completed")

    def _process_clusters(self,cluster_data,barcode_index_map, expr_matrix, all_barcodes):
        for cluster, cluster_barcodes in tqdm(cluster_data.items(),desc="Processing clusters"):
            self.logger.info(f"Processing cluster: {cluster} with {len(cluster_barcodes)} barcodes")
            
            sub_expr, sub_barcodes = self.matrix_builder.extract_cluster_expression(barcode_index_map, expr_matrix, all_barcodes, cluster_barcodes)

            if sub_expr is None:
                self.logger.warning(f"Cluster {cluster}: No valid barcodes.")
                continue

            # save cluster expression
            self.gene_writer.save_cluster_expression(cluster, sub_expr, sub_barcodes)


