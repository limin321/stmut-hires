import os
import logging

import pandas as pd


ANNOTATE_REQUIRED_COLS = {"cluster", "annotate"}
BULK_CNV_REQUIRED_COLS = {"arms", "gainloss"}


class InputValidator:
    """
    Validate pipeline runtime state, inputs, metadata,
    and required outputs for downstream execution.
    """

    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)

    # ==========================================================
    # Public APIs
    # ==========================================================
    def validate_run(self):
        """Validate full pipeline execution requirements."""
        self.logger.info(
            "Validating full pipeline inputs..."            
        )
        self._ensure_output_dir()
        if not self.config.dry_run:
            self._validate_canonical_inputs()
            
        self._validate_pipeline_parameters()
        self._validate_cnv_inputs()
        self.logger.info(
            "Full pipeline validation completed successfully."
        )


    def validate_call_cnv(self):
        """Validate Step 6 standalone execution requirements."""
        self._validate_metadata()
        self._validate_previous_pipeline_outputs()
        self._validate_cnv_inputs()
        self._validate_canonical_inputs()
        self.logger.info(
            "call-cnv validation completed successfully."
        )

    # ==========================================================
    # Output directory
    # ==========================================================
    def _ensure_output_dir(self):
        if self.config.dry_run:
            self.logger.info(
                f"[DRY-RUN] Would create output dir: "
                f"{self.config.output_dir}"
            )

            return

        os.makedirs(
            self.config.output_dir, 
            exist_ok=True
        )

        self.logger.info(
            f"Output directory ready: "
            f"{self.config.output_dir}"
        )

    # ==========================================================
    # Metadata validation
    # ==========================================================
    def _validate_metadata(self):
        """ 
        Ensure metadata.json exists for downstream execution.
        """    
        if not self.config.metadata_manager:
            raise ValueError(
                "metadata_manager is not configured."
            )
        metadata_path = (
            self.config.metadata_manager.metadata_path
        )

        if not os.path.exists(metadata_path):
            raise FileNotFoundError(
                f"Pipeline metadata file not found."
                f"{metadata_path}"
            )

        self.logger.info(
            f"Metadata validated: {metadata_path}"
        )

    # ==========================================================
    # Canonical pipeline inputs
    # ==========================================================
    def _validate_canonical_inputs(self):
        """ 
        Validate canonical runtime inputs resolved from metadata.json
        """  
        self.logger.info(
            "Validating canonical pipeline inputs ..."
        )

        self._validate_required_file(
            filepath=self.config.clusterf,
            arg_name="cluster_file"
        )

        self._validate_required_file(
            filepath=self.config.exp_h5,
            arg_name="exp_h5"
        )

        self._validate_required_file(
            filepath=self.config.spatial_file,
            arg_name="spatial_file"
        )

        self._validate_required_file(filepath=self.config.bed_file, arg_name="bed_file")
        self._validate_required_file(filepath=self.config.gene_info, arg_name="gene_info")
        self._validate_required_file(filepath=self.config.corr_file, arg_name="corr_file")
        self.logger.info(
            "Canonical inputs validated successfully."
        )

    # ==========================================================
    # Pipeline runtime parameters
    # ==========================================================
    def _validate_pipeline_parameters(self):

        self.logger.info(
            "Validating runtime parameters..."
        )

        self._validate_type(
            self.config.cutoff,
            int,
            "cutoff"
        )

        self._validate_type(
            self.config.window,
            int,
            "window"
        )

        self._validate_type(
            self.config.cores,
            int,
            "cores"
        )

        if self.config.num_processes is not None:

            self._validate_type(
                self.config.num_processes,
                int,
                "num_processes"
            )

        self.logger.info(
            "Runtime parameters validated successfully."
        )    

    # ==========================================================
    # CNV Step inputs
    # ==========================================================
    def _validate_cnv_inputs(self):

        self.logger.info(
            "Validating CNV step inputs..."
        )

        # ------------------------------------------------------
        # annotate_file (required)
        # ------------------------------------------------------
        self._validate_required_file(
            filepath=self.config.annotate_csv,
            arg_name="annotate_file"
        )

        self._check_csv_columns(
            filepath=self.config.annotate_csv,
            required_cols=ANNOTATE_REQUIRED_COLS,
            arg_name="annotate_file"
        )

        self.logger.info(
            f"Annotation file validated: "
            f"{self.config.annotate_csv}"
        )

        # ------------------------------------------------------
        # bulkCNV_file (optional)
        # ------------------------------------------------------
        if self.config.bulk_csv:

            self._validate_required_file(
                filepath=self.config.bulk_csv,
                arg_name="bulkCNV_file"
            )

            self._check_csv_columns(
                filepath=self.config.bulk_csv,
                required_cols=BULK_CNV_REQUIRED_COLS,
                arg_name="bulkCNV_file"
            )

            self.logger.info(
                f"Bulk CNV file validated: "
                f"{self.config.bulk_csv}"
            )

        else:

            self.logger.info(
                "No bulk CNV file provided."
            )

    # ==========================================================
    # Previous pipeline outputs
    # ==========================================================

    def _validate_previous_pipeline_outputs(self):
        """
        Validate outputs required for re-running Step 6.
        """

        wtcnr_dir = os.path.join(
            self.config.output_dir,
            "wtcnr"
        )

        if not os.path.exists(wtcnr_dir):

            raise FileNotFoundError(
                f"Required wtcnr directory not found: "
                f"{wtcnr_dir}. "
                f"Please run the full pipeline first."
            )

        cnr_files = [
            f for f in os.listdir(wtcnr_dir)
            if f.endswith(".cnr")
        ]

        if not cnr_files:

            raise FileNotFoundError(
                f"No .cnr files found in: "
                f"{wtcnr_dir}"
            )

        self.logger.info(
            f"Validated {len(cnr_files)} "
            f"weighted CNR files."
        )

    # ==========================================================
    # Generic reusable validators
    # ==========================================================

    def _validate_required_file(
        self,
        filepath,
        arg_name
    ):

        if not filepath:

            raise ValueError(
                f"{arg_name} is required."
            )

        if not os.path.exists(filepath):

            raise FileNotFoundError(
                f"{arg_name} not found: {filepath}"
            )

    def _validate_type(
        self,
        value,
        expected_type,
        param_name
    ):

        if not isinstance(value, expected_type):

            raise TypeError(
                f"{param_name} must be "
                f"{expected_type}, got {type(value)}"
            )

    def _check_csv_columns(
        self,
        filepath,
        required_cols,
        arg_name
    ):
        """
        Read only CSV header and validate columns.
        """

        actual_cols = set(
            pd.read_csv(filepath, nrows=0).columns
        )

        missing = required_cols - actual_cols

        if missing:

            raise ValueError(
                f"{arg_name} missing required columns: "
                f"{missing}. "
                f"Expected: {required_cols}. "
                f"Observed: {actual_cols}"
            )








