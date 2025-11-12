"""Gene Data Writer """
import os
import pandas as pd
import logging

class GeneDataWriter:
    """Handles writing gene data and expression file """
    def __init__(self, output_dir, dry_run=False):
        self.output_dir = os.path.join(output_dir,"cluster_exp")
        self.dry_run = dry_run
        self.logger = logging.getLogger(__name__)

        if not self.dry_run:
            os.makedirs(self.output_dir, exist_ok=True)

    def save_gene_ids(self, genes):
        """Save gene ensembl ID"""
        output_path = os.path.join(self.output_dir, "ensembl.csv")
        
        if self.dry_run:
            self.logger.info(f"[DRY-RUN] Would save {len(genes)} gene IDs to: {output_path}")
        else:
            pd.DataFrame({'gene_id': genes}).to_csv(output_path, index=False)
            self.logger.info(f"Saved {len(genes)} gene IDs to: {output_path}")

    def save_cluster_expression(self, cluster_name, cluster_exp, cluster_barcodes):
        """Save cluster gene expression matrix to parquet """
        output_path = os.path.join(self.output_dir, f"{cluster_name}.parquet")
        
        if self.dry_run:
            self.logger.info(f"[DRY-RUN] Would save cluster '{cluster_name}' expression data to: {output_path}")
            self.logger.info(f"[DRY-RUN] Cluster '{cluster_name}' shape: {cluster_exp.shape}")
        else:
            df = pd.DataFrame(cluster_exp, columns=cluster_barcodes)
            df.to_parquet(
                output_path,
                index=False,
                engine='pyarrow',
                compression='snappy'
            )
            self.logger.info(f"Saved cluster '{cluster_name}' expression data to: {output_path}")
