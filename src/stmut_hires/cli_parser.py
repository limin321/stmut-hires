import argparse

class CommandLineParser:
    """Responsible for parsing command line arguments"""

    @staticmethod
    def parse_args():
        parser = argparse.ArgumentParser(description="CNV Analysis Pipeline.")
        subparsers = parser.add_subparsers(dest="command", help="Available commands")

        # 1. Global Parent: Required for All commands (Step 1-6)
        global_parent = argparse.ArgumentParser(add_help=False)
        global_parent.add_argument("--output_dir", required=True, help="Output directory.")
        global_parent.add_argument("--cluster_file", required=True, help="Barcodes cluster info.(csv);Graph-based.csv downloaded from Loupe Browser.")
        global_parent.add_argument("--dry_run", action="store_true", help="Dry run mode")    
        
        # 2. Step 1-5 Parent: Required only for Full pipeline
        step1_5_parent = argparse.ArgumentParser(add_help=False)
        step1_5_parent.add_argument("--exp_h5", required=True, help="Gene expression matrix (h5)")
        step1_5_parent.add_argument("--spatial_file", required=True, help="Barcodes spatial coordinates.(./spatial/tissue_positions.parquet)")
        step1_5_parent.add_argument("--manual_cutoff", required=False, default= None, type=int,
                            help="Cutoff to filter-out barcodes with low gene counts (INT, default: %(default)s). You either set this parameter or `--bw_method` to filter-out low-quality barcodes.")
        step1_5_parent.add_argument("--num_processes", required=False, default=None,type=int, 
                            help="The number of processes used for parallel (default: %(default)s)")
        step1_5_parent.add_argument("--cutoff", required=False, default=1000,type=int, 
                            help="The min number of genes in a new spot as grouping cutoff. (default: %(default)s)")
        step1_5_parent.add_argument("--window", required=False, default=100,type=int, 
                            help="The nearest-neighbor spots for selecting grouping candidates. (default): %(default)s")
        step1_5_parent.add_argument("--cores", default=4, type=int, help="Number of cores to run weighted-median parallelly")
        step1_5_parent.add_argument("--bw_method", required=False, default=0.1, type=float,
                            help="Bandwidth scaling factor for KDE valley detection (e.g., 0.1 to 0.9). For it to work, `--manual_cutoff` needs to be default value None.")
        
        # 3. Step 6 Parent: Required for both 'run' and 'call-cnv'
        step6_parent = argparse.ArgumentParser(add_help=False)
        step6_parent.add_argument("--annotate_file", required=False, default=None, type=str, help="Two-column clusters annotated csv file.")
        step6_parent.add_argument("--bulkCNV_file", required=False, default=None, type=str, help="[Optional] Two-column csv storing bulk-CNV info.")
        step6_parent.add_argument("--pmtimes", default=5, type=int, help="The number of permutation times.")
        step6_parent.add_argument("--ncluster", default=6, type=int, help="Number of clusters for CNV plot.")
        step6_parent.add_argument("--distance_metric", type=str, help="Distance metric for clustering")
        step6_parent.add_argument("--linkage_method", type=str, help="Linkage method for hierarchical clustering")
        
        # Sub-command: run (Full pipeline)
        subparsers.add_parser(
            "run", 
            parents=[global_parent,step1_5_parent,step6_parent], 
            help="Run full pipeline."
        )
        
        # Sub-command: call-cnv (Step 6)
        subparsers.add_parser(
            "call-cnv", 
            parents=[global_parent,step6_parent], 
            help="Run Step 6: calling CNV..only")

        # Independent clean sub-command
        clean_step = argparse.ArgumentParser(add_help=False)
        clean_step.add_argument("--output_dir", required=True, help="Output directory to clean.")
        
        subparsers.add_parser(
            "clean",
            parents=[clean_step],
            help="Clean intermediate output folders, keeping only figures/ and tables/.",

        )
        
        return parser.parse_args()


