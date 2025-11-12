#! /usr/bin/env python3

import os
import sys
import logging
import glob
import pyarrow
import pyarrow.parquet as pq

# Add the workspace root to Python path (two levels up from tools/)
workspace_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, workspace_root)

# Multiprocessing
from parallel.cli_parser import CommandLineParser
from config import MergerConfig
from split_cluster_expr.cluster_processor import ClusterExpressionProcessor
from parallel.parallel_processor import ParallelClusterProcessor


def setup_logging():
    """Setup basic logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def verify_parquet_files(file_pattern, logger):
    """
    Quickly checks if all parquet files structure are valid by attempting to read their schema.
    Raise an IOError if any file is corrupted or none are found.
    """
    files = glob.glob(file_pattern)
    if not files:
        dir_path = os.path.dirname(file_pattern)
        raise IOError(f"NO parquet files are detected in directory: ${dir_path}")
    
    logger.info(f"Verifying {len(files)} parquet files...")
    for file_path in files:
        try:
            pq.read_schema(file_path)
        except(IOError, pyarrow.lib.ArrowInvalid) as e:
            logger.error(f"Failed verification for file: {file_path}")
            raise IOError(f"Corrupted Parquet file detected: {file_path}. Error: {e}")
    logger.info(f"All {len(files)} parquet files verified successfully.")
    return True


def main():    
    args = CommandLineParser.parse_args()
    # Setup logging
    logger = setup_logging()

    try:
        config = MergerConfig(
            clusterf=args.cluster_file,
            exp_h5=args.exp_h5,
            output_dir=args.output_dir,
            dry_run=args.dry_run,
            spatial_file=args.spatial_file,
            cutoff=args.cutoff,
            window=args.window,
            num_processes=args.num_processes
        )
    except (FileNotFoundError, TypeError, NotADirectoryError) as e:
        logger.error(f"Configuration Error: {e}")
        sys.exit(1)


    # Define the directory where intermediate cluster parquet files are stored.
    intermediate_cluster_dir = os.path.join(config.output_dir, "cluster_exp")
    expression_file_pattern = os.path.join(intermediate_cluster_dir, "Cluster*.parquet")
    ensembl_path = os.path.join(intermediate_cluster_dir, "ensembl.csv")

    files_exist = glob.glob(expression_file_pattern)  

    # Set the config paths right away to where the files *should* be
    config.set_step1_outputs(
        expression_file_path = expression_file_pattern,
        ensembl_file_path = ensembl_path
    )

    if files_exist:
        logger.info(f"Intermediate cluster files found in '{intermediate_cluster_dir}'.")

        try:
            verify_parquet_files(expression_file_pattern, logger)

        except IOError as e:
            # This handles both "No files found" and "Corrupted files detected" errors
            logger.error(e)
            logger.info("Corrupt files detected. Forcing re-run of Step 1 processing.")
            # Optional: Clean up existing broken files before the re-run starts
            for f in glob.glob(expression_file_pattern) + [ensembl_path]:
                 if os.path.exists(f):
                     os.remove(f)
                     logger.info(f"Deleted broken file: {f}")

            files_exist = False

    if not files_exist:
        logger.info("Starting step 1 processing...")
        processor = ClusterExpressionProcessor(config, dry_run=args.dry_run)
        processor.process()

    # The crucial change: The validation should happen *here*, 
    # after we are certain that either existing files passed verification, 
    # or new files have just been created by processor.process().
    
    # Then call the config obj defined above, it inherites other classes in the config.py file.
    config.validate_inputs()

    if args.dry_run:
        logger.info("Starting dry-run mode - no files will be written")
    else:
        logger.info("Starting normal processing mode")
    
    logger.info("Starting Step 2 (ParallelClusterProcessor/Merger)...")
    # Process clusters in parallel
    parallel_processor = ParallelClusterProcessor(config=config)
    parallel_processor.parallel_cluster_merger()
    
if __name__ == "__main__":
    main()

