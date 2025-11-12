import argparse

class CommandLineParser:
    """Responsible for parsing command line arguments"""

    @staticmethod
    def parse_args():
        parser = argparse.ArgumentParser(description="Grouping barcodes within each cluster.")
        parser.add_argument("--exp_h5", required=True, help="Gene expression matrix (h5)")
        parser.add_argument("--cluster_file", required=True, help="Barcodes cluster info.(csv)")
        parser.add_argument("--spatial_file", required=True, help="Barcodes spatial coordinates.(./spatial/tissue_positions.parquet)")
        parser.add_argument("--output_dir", required=True, help="Output directory.")
        parser.add_argument("--num_processes", required=False, default=None,type=int, 
                            help="The number of processes used for parallel (default: %(default)s)")
        parser.add_argument("--cutoff", required=False, default=1000,type=int, 
                            help="The min number of genes in a new spot as grouping cutoff. (default: %(default)s)")
        parser.add_argument("--window", required=False, default=100,type=int, 
                            help="The nearest-neighbor spots for selecting grouping candidates. (default): %(default)s")
        parser.add_argument("--dry_run", action="store_true", help="Dry run mode")    
        return parser.parse_args()


