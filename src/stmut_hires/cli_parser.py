import argparse

class CommandLineParser:
    """Parse and validate command line arguments"""

    SUPPORTED_PLATFORMS = [
        "visiumhd",
        "atera"
    ]

    @staticmethod
    def parse_args():
        parser = argparse.ArgumentParser(description="CNV Analysis Pipeline.")
        subparsers = parser.add_subparsers(
            dest="command", 
            required = True,
            help="Available commands"
        )

        # =========================================
        # Shared arguments
        # =========================================
        # 1. Global Parent: Required for All commands (Step 1-6)
        io_parent = argparse.ArgumentParser(add_help=False)
        io_parent.add_argument("--output_dir", required=True, help="Pipeline output directory.")
        io_parent.add_argument("--dry_run", action="store_true", help="Dry run mode. Validate inputs only.")


        platform_parent = argparse.ArgumentParser(add_help=False)
        platform_parent.add_argument("--input_mode", required=True, choices=CommandLineParser.SUPPORTED_PLATFORMS, help="Input platform type")
        platform_parent.add_argument("--exp_h5",help="Spatial gene expression matrix h5")

        # =======================================================
        # Canonical pipeline inputs (used directly by pipeline)
        # Used VisiumHD outs as standard formats.
        # =======================================================
        canonical_parent = argparse.ArgumentParser(add_help=False)
        canonical_parent.add_argument("--spatial_file",help="[VisiumHD] Barcodes spatial coordinates.(./spatial/tissue_positions.parquet)")
        canonical_parent.add_argument("--cluster_file", help="[VisiumHD] Barcodes cluster info.(csv);Graph-based.csv downloaded from Loupe Browser.")

        # =======================================================
        # Atera raw inputs 
        # =======================================================
        atera_parent = argparse.ArgumentParser(add_help=False)
        atera_parent.add_argument("--cells_parquet", help="[Atera] cells.parquet")
        atera_parent.add_argument("--analysis_zarr_zip", help="[Atera] analysis.zarr.zip")

        # =======================================================
        # Pipeline parameters
        # =======================================================
        pipeline_parent = argparse.ArgumentParser(add_help=False)
        pipeline_parent.add_argument("--manual_cutoff", required=False, default= None, type=int,
                            help="Cutoff to filter-out barcodes with low gene counts (INT, default: %(default)s). You either set this parameter manually or provide `--bw_method` to automatically predict a value.")

        pipeline_parent.add_argument("--num_processes", required=False, default=None,type=int, 
                            help="The number of processes used for parallel (default: %(default)s)")
        pipeline_parent.add_argument("--cutoff", required=False, default=1000,type=int, 
                            help="The min number of genes in a new spot as grouping cutoff. (default: %(default)s)")
        pipeline_parent.add_argument("--window", required=False, default=100,type=int, 
                            help="The nearest-neighbor spots for selecting grouping candidates. (default: %(default)s)")
        pipeline_parent.add_argument("--cores", default=4, type=int, help="Number of cores to run weighted-median parallelly")
        pipeline_parent.add_argument("--bw_method", required=False, default=0.1, type=float,
                            help="Bandwidth scaling factor for KDE valley detection (e.g., 0.1 to 0.9). For it to work, `--manual_cutoff` needs to be default value None.")
        pipeline_parent.add_argument("--smooth_method", required=False, default="arm", choices=["arm", "local"],
                    help="Weighted-median smoothing method: 'arm' or 'local' (default: %(default)s)")
        pipeline_parent.add_argument("--target_weight", required=False, default=25, type=int,
                            help="Target weight for local smoothing (only used when --smooth_method=local, default: %(default)s)")

        # =======================================================
        # CNV step parameters
        # =======================================================
        cnv_parent = argparse.ArgumentParser(add_help=False)
        cnv_parent.add_argument("--annotate_file", required=True, default=None, type=str, help="Two-column clusters annotated csv file.")
        cnv_parent.add_argument("--bulkCNV_file", required=False, default=None, type=str, help="[Optional] Two-column csv storing bulk-CNV info.")
        cnv_parent.add_argument("--pmtimes", default=5, type=int, help="[Optional] The number of permutation times.")
        cnv_parent.add_argument("--ncluster", default=6, type=int, help="Number of clusters for CNV plot, (default: %(default)s).")
        cnv_parent.add_argument("--distance_metric",default="euclidean" ,type=str, help="Distance metric for clustering")
        cnv_parent.add_argument("--linkage_method", type=str,default="ward", help="Linkage method for hierarchical clustering")

        # =======================================================
        # RUN command
        # =======================================================
        # Run full pipeline
        subparsers.add_parser(
            "run",
            parents=[
                io_parent,
                platform_parent,
                canonical_parent,
                atera_parent,
                pipeline_parent,
                cnv_parent
            ],
            help="Run full pipeline"
        )

        
        # Sub-command: call-cnv (Step 6)
        subparsers.add_parser(
            "call-cnv", 
            parents=[io_parent, cnv_parent], 
            help="Run Step 6: calling CNV..only")

        # Independent clean sub-command
        clean_parser = subparsers.add_parser("clean", help="Clean intermediate files")
        clean_parser.add_argument("--output_dir", required=True, help="Output directory to clean.")
        
        args = parser.parse_args()
        CommandLineParser.validate_args(args)
        return args

    @staticmethod
    def validate_args(args):

        if args.command != "run":
            return
    
        # VisiumHD mode
        if args.input_mode == "visiumhd":
            required = {
                "--exp_h5": args.exp_h5,
                "--spatial_file": args.spatial_file,
                "--cluster_file": args.cluster_file
            }

            missing = [
                k for k, v in required.items()
                if v is None
            ]

            if missing:
                raise ValueError(
                    f"VisiumHD mode missing required args: {missing}"
                )

        # Atera mode
        elif args.input_mode == "atera":
            required = {
                "--exp_h5": args.exp_h5,
                "--analysis_zarr_zip": args.analysis_zarr_zip,
                "--cells_parquet": args.cells_parquet
            }

            missing = [
                k for k, v in required.items()
                if v is None
            ]

            if missing:
                raise ValueError(
                    f"Atera mode missing required args: {missing}"
                )

            

        


