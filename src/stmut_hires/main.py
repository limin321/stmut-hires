#! /usr/bin/env python3

import os
import sys
import logging
import glob
import subprocess

import pyarrow
import pyarrow.parquet as pq

# # Add the workspace root to Python path (two levels up from tools/)
# workspace_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# sys.path.insert(0, workspace_root)

# Multiprocessing
from stmut_hires.cli_parser import CommandLineParser
from stmut_hires.config import CNVAnalysisConfig
from stmut_hires.data_io.validators import InputValidator
from stmut_hires.split_cluster_expr.cluster_expr_generator import ClusterExpressionProcessor
from stmut_hires.parallel_cluster_runner.parallel_processor import ParallelClusterProcessor
from stmut_hires.weighted_median.parallel_wtmedian import WorkflowOrchestrator
from stmut_hires.call_cnv.callCNV_workflow import CNVCallPlotWorkflow
from stmut_hires.parallel_cluster_runner.output_manager import OutputDirManager
from stmut_hires.data_io.output_cleaner import OutputCleaner
from stmut_hires.adapters.atera_adapter import AteraAdapter
from stmut_hires.metadata.metadata import PipelineMetadataManager




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
        raise IOError(f"NO parquet files are detected in directory: {dir_path}")
    
    logger.info(f"Verifying {len(files)} parquet files...")

    for file_path in files:
        try:
            pq.read_schema(file_path)
        except(IOError, pyarrow.lib.ArrowInvalid) as e:
            logger.error(f"Failed verification for file: {file_path}")
            raise IOError(f"Corrupted Parquet file detected: {file_path}. Error: {e}")
    logger.info(f"All {len(files)} parquet files verified successfully.")
    return True

"""  
Design Structure:
CLI
  ↓
platform adapter
  ↓
metadata.json
  ↓
config.resolve_canonical_inputs()
  ↓
pipeline



cli_parser/
    only parses args

adapters/
    converts platform formats

metadata/
    stores pipeline state

config/
    stores runtime resolved configuration

validators/
    validates resolved runtime state

workflows/
    execute business logic


"""

def prepare_inputs_for_run(args_dict, logger):
    """ 
    Resolve canonical + raw inputs. Does NOT write metadata.json.
    """

    input_mode = args_dict["input_mode"]

    # Prepare canonical inputs for the pipeline regardless of platform
    if input_mode == "visiumhd":
        canonical_inputs = {
            "cluster_file": args_dict["cluster_file"],
            "spatial_file": args_dict["spatial_file"],
            "exp_h5": args_dict["exp_h5"]
        }

        raw_inputs = canonical_inputs.copy()

    elif input_mode == "atera":
        logger.info("Running Atera adapter conversion...")
        atera_converter = AteraAdapter(args_dict["output_dir"], dry_run=args_dict.get("dry_run", False))
        cluster_file = atera_converter.get_graph_based(
            args_dict["analysis_zarr_zip"], args_dict["cells_parquet"]
        )
        spatial_file = atera_converter.get_tissue_position_parquet(
            args_dict["cells_parquet"]
        )

        canonical_inputs = {
            "cluster_file": cluster_file,
            "spatial_file": spatial_file,
            "exp_h5": args_dict["exp_h5"],
        }
        raw_inputs = {
            "cells_parquet": args_dict["cells_parquet"],
            "analysis_zarr_zip": args_dict["analysis_zarr_zip"],
            "exp_h5": args_dict["exp_h5"],
        }

    
    else:
        raise ValueError(f"Unsupported input mode: {input_mode}")

    # Create metadata manager
    metadata_manager = PipelineMetadataManager(
        output_dir = args_dict.get("output_dir")
    )
    return metadata_manager, canonical_inputs, raw_inputs


def build_config(args_dict, metadata_manager, canonical_inputs):
    """ 
    Build runtime config object.
    Canonical inputs are resolved from metadata.json
    """
    # Use vars(args).get() to safely handle arguments that might be missing
    # depending on which subcommand (run vs call-cnv) was used.
    config = CNVAnalysisConfig(
        canonical_inputs=canonical_inputs,
        metadata_manager=metadata_manager,

        filter_cutoff = args_dict.get('manual_cutoff'),
        bw_method = args_dict.get('bw_method'),

        num_processes=args_dict.get('num_processes'),
        cutoff=args_dict.get('cutoff', 1000),
        window=args_dict.get('window', 100),

        output_dir=args_dict.get('output_dir'),
        dry_run=args_dict.get('dry_run'),
        cores=args_dict.get('cores', 4),

        annotate_csv=args_dict.get('annotate_file'),
        bulk_csv=args_dict.get('bulkCNV_file'),

        ncluster=args_dict.get('ncluster', 6),
        pmtimes=args_dict.get('pmtimes', 5),
        distance_metric=args_dict.get(
            'distance_metric', 
            'euclidean'
        ),
        linkage_method=args_dict.get(
            'linkage_method', 
            'ward'
        ),
    )

    return config

def main():    
    args = CommandLineParser.parse_args()
    logger = setup_logging()

    # Handle clean before config is built — it doesn't need building config for clean step only.
    if args.command == "clean":
        OutputCleaner(args.output_dir).clean()
        return

    args_dict = vars(args)

    # ==========================================================
    # Resolve inputs + build config
    # (metadata.json is NOT written yet for `run` — done after validation)
    # ==========================================================
    try:
        if args.command == "run":
            metadata_manager, canonical_inputs, raw_inputs = prepare_inputs_for_run(
                args_dict, logger
            )

        elif args.command == "call-cnv":
            metadata_manager = PipelineMetadataManager(
                output_dir=args_dict["output_dir"]
            )
            md = metadata_manager.load_metadata()
            # Re-running Step 6: pull canonical inputs from the existing metadata.json
            canonical_inputs = md["canonical_inputs"]
            for k, v in md.get("runtime_params", {}).items():
                if args_dict.get(k) is None:
                    args_dict[k] = v
            raw_inputs = None  # unused for call-cnv

        else:
            raise ValueError(f"Unsupported command: {args.command}")

        config = build_config(args_dict, metadata_manager, canonical_inputs)

    except (
        FileNotFoundError,
        TypeError,
        NotADirectoryError,
        ValueError,
        KeyError,
    ) as e:
        logger.error(f"Configuration Error: {e}")
        sys.exit(1)

    # ==========================================================
    # Run Pipeline
    # ==========================================================
    if args.command == "run":
        logger.info("Running full pipeline ...")

        # Define the directory where intermediate cluster parquet files are stored.
        intermediate_cluster_dir = os.path.join(config.output_dir, "cluster_exp")
        expression_file_pattern = os.path.join(intermediate_cluster_dir, "Cluster*.parquet")
        ensembl_path = os.path.join(intermediate_cluster_dir, "ensembl.csv")

        config.set_step1_outputs(
            expression_file_path=expression_file_pattern,
            ensembl_file_path=ensembl_path,
        )

        # Validate all raw inputs before doing any work.
        InputValidator(config).validate_run()

        if config.dry_run:
            logger.info("Dry-run complete. Inputs validated successfully.")
            return

        # Only persist metadata after inputs are confirmed good.
        metadata_manager.initialize_metadata(
            platform=args_dict["input_mode"],
            canonical_inputs=canonical_inputs,
            raw_inputs=raw_inputs,
            runtime_params={
                "bw_method": config.bw_method,
                "manual_cutoff": config.filter_cutoff,
                "cutoff": config.cutoff,
                "window": config.window,
                "distance_metric": config.distance_metric,
                "linkage_method": config.linkage_method,
                "ncluster": config.ncluster,
                "pmtimes": config.pmtimes,
            },
            completed_steps=["adapter_conversion"],
        )
        logger.info("Metadata initialized successfully.")

        # STEP 1 ...
        files_exist = glob.glob(expression_file_pattern)

        if files_exist:
            logger.info(f"Intermediate cluster files found in '{intermediate_cluster_dir}'.")

            try:
                verify_parquet_files(expression_file_pattern, logger)

            except IOError as e:
                logger.error(e)
                logger.info(
                    "Corrupt files detected. "
                    "Re-running of Step 1."
                )

                for f in glob.glob(expression_file_pattern) + [ensembl_path]:
                    if os.path.exists(f):
                        os.remove(f)
                        logger.info(f"Deleted broken file: {f}")

                files_exist = False

        if not files_exist:
            logger.info("Starting step 1 processing...")
            processor = ClusterExpressionProcessor(config)
            processor.process()
        metadata_manager.update_completed_steps("step1")

        # Step 2 ...
        logger.info("Starting Step 2 (ParallelClusterProcessor/Merger)...")

        output_creator = OutputDirManager(config.output_dir)
        txt_output_dir, cnr_dir, wtcnr_dir, cdt_dir,_,_ = output_creator.create_output_dir()

        summary_exist = glob.glob(os.path.join(
            config.output_dir, 
            "cluster_summary",
            "*barcode_grouping_info.csv"
        ))

        if not summary_exist:
            # Process clusters in parallel
            parallel_processor = ParallelClusterProcessor(config=config)
            parallel_processor.parallel_cluster_merger()
            output_creator.move_csvs_to_summary() # move merger summary to cluster_summary dir
        else:
            logger.info("Step 2 already completed. ")

        metadata_manager.update_completed_steps("step2")

        # Step 3 ...
        logger.info("Starting Step 3:  "
                "Generate cnr file for each barcode ...")

        cnr_exist = glob.glob(os.path.join(cnr_dir, "*.cnr"))

        if cnr_exist:
            logger.info("cnr files are generated.")
            logger.warning(f"If you want to regenerate cnr files, make sure {cnr_dir} is empty.")
        
        if not cnr_exist:
            patch_code = (
                "import pandas as pd; "
                "pd.DataFrame.iteritems = pd.DataFrame.items; "
                "pd.Series.iteritems = pd.Series.items; "
                "import sys; "
                "from cnvlib.cnvkit import main; "
                "sys.exit(main())"
            )

            txt_files = glob.glob(os.path.join(txt_output_dir, "*.txt"))

            logger.info(f"Running import-rna on {len(txt_files)} files...")

            full_cmd = (
                f"ulimit -s unlimited && cd {txt_output_dir} && "
                f'python -c "{patch_code}" import-rna '
                f"-f counts -g {config.gene_info} -c {config.corr_file} "
                f"--output-dir {cnr_dir} " 
                f"-o {cnr_dir}/output.txt "
                f"$(ls *.txt)"  # Command substitution to pass all files at once
            )
            subprocess.run(
                full_cmd, 
                shell=True, 
                check=True
            )
        else:
            logger.info("cnr files already exist.")

        metadata_manager.update_completed_steps("step3")

        # Step 4 ...
        logger.info("Starting step 4: weighted_median calculation...")

        wtcnr_exist = glob.glob(os.path.join(wtcnr_dir, "*.cnr"))

        if not wtcnr_exist:
            orchestrator = WorkflowOrchestrator(config=config)
            orchestrator.run_parallel_workflow()
        else:
            logger.info("Weighted cnr files already exist.")

        metadata_manager.update_completed_steps("step4")

        # Step 5 ...
        logger.info("Starting step 5: export cdt file...")

        cnr_files_list = glob.glob(os.path.join(wtcnr_dir, "*.cnr"))

        export_shell_cmd = (
            f"ulimit -s unlimited && "
            f"cd {wtcnr_dir} && "
            f"cnvkit.py export cdt *.cnr -o {cdt_dir}/grpWt.cdt"
        )

        logger.info(f"Exporting CDT for {len(cnr_files_list)} files...")
        
        try:
            subprocess.run(
                export_shell_cmd, 
                shell=True, 
                check=True
            )
            logger.info("Export successful.")

        except subprocess.CalledProcessError as e:
            logger.error(f"CDT export failed: {e}")

        metadata_manager.update_completed_steps("step5")

        # Step 6 ...
        logger.info("Starting step 6: calling CNV...")

        workflow = CNVCallPlotWorkflow(config=config)
        workflow.cnv_workflow(
            ncluster=config.ncluster,
            bulk_csv=config.bulk_csv,
            pmtimes=config.pmtimes,
            distance_metric=config.distance_metric,
            linkage_method=config.linkage_method,
            )

        metadata_manager.update_completed_steps("step6")
# =======================
# call-cnv only
# =======================
    elif args.command == "call-cnv":
        logger.info("Re-running Step 6: calling CNV ...")
        InputValidator(config).validate_call_cnv()
        workflow = CNVCallPlotWorkflow(config=config)
        workflow.cnv_workflow(
            ncluster=config.ncluster,
            bulk_csv=config.bulk_csv,
            pmtimes=config.pmtimes,
            distance_metric=config.distance_metric,
            linkage_method=config.linkage_method,
        ) 
        metadata_manager.update_completed_steps("step6")

if __name__ == "__main__":
    main()


""" 
indir="/stomics_data/liminData/Visium/stmut_python/BD17_bin8_inputs"
outdir="/stomics_data/liminData/Visium/stmut_python/outs1"

# step1&2
python ./src/stmut_hires/main.py \
    --exp_h5 ${indir}/filtered_feature_bc_matrix.h5 \
    --cluster_file ${indir}/Graph-Based.csv \
    --spatial_file ${indir}/spatial/tissue_positions.parquet \
    --output_dir ${outdir} \
    --cutoff 3000 \
    --num_processes 3 \
    --annotate_file ${indir}/annotate.csv \
    --bulkCNV_file ${indir}/bulkCNV.csv \
    --cores 6 \
    --pmtimes 5 \
    --ncluster 6

# set up cnvkit
git clone https://github.com/etal/cnvkit
cd cnvkit/
pip install -e .
cnvkit.py import-rna --help

# to add sub_command, modify config, cli_parser, main three files.
python ../src/stmut_hires/main.py call-cnv \
    --cluster_file ${indir}/Graph-Based.csv \
    --output_dir ${outdir} \
    --annotate_file ${indir}/annotate.csv \
    --bulkCNV_file ${indir}/bulkCNV.csv \
    --pmtimes 5 \
    --ncluster 6
"""


