import os
import logging

import pandas as pd


ANNOTATE_REQUIRED_COLS = {"cluster", "annotate"}
BULK_CNV_REQUIRED_COLS = {"arms", "gainloss"}


class InputValidator:
    """Validates presence and correctness of all pipeline input files."""

    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)

    def validate_run(self):
        """Validate all raw input files for the full 'run' pipeline."""
        self._ensure_output_dir()
        self._validate_step1()
        self._validate_step2()
        self._validate_step6()

    def validate_call_cnv(self):
        """Validate inputs for the 'call-cnv' subcommand."""
        self._ensure_output_dir()
        self._validate_step6()

    def _ensure_output_dir(self):
        if not self.config.dry_run:
            os.makedirs(self.config.output_dir, exist_ok=True)
            self.logger.info(f"Created output dir: {self.config.output_dir}")
        else:
            self.logger.info(f"[DRY-RUN] Would create output dir: {self.config.output_dir}")

    def _validate_step1(self):
        self.logger.info("Validating step 1 input files...")
        if not self.config.clusterf:
            raise ValueError("--cluster_file is required.")
        if not os.path.exists(self.config.clusterf):
            raise FileNotFoundError(f"Cluster file NOT FOUND: {self.config.clusterf}")

        if not self.config.exp_h5:
            raise ValueError("--exp_h5 is required.")
        if not os.path.exists(self.config.exp_h5):
            raise FileNotFoundError(f"H5 file NOT FOUND: {self.config.exp_h5}")

        self.logger.info("Step 1 inputs are valid.")

    def _validate_step2(self):
        self.logger.info("Validating Step 2 specific inputs...")
        if not self.config.spatial_file:
            raise ValueError("--spatial_file is required.")
        if not os.path.exists(self.config.spatial_file):
            raise FileNotFoundError(f"Spatial input file not found: {self.config.spatial_file}")
        if not isinstance(self.config.cutoff, int):
            raise TypeError(f"cutoff must be an integer, got: {type(self.config.cutoff)}")
        if not isinstance(self.config.window, int):
            raise TypeError(f"window must be an integer, got: {type(self.config.window)}")
        self.logger.info("Step 2 inputs are valid.")

    def _validate_step6(self):
        # annotate.csv is required
        if not self.config.annotate_csv:
            raise ValueError("--annotate_file is required.")
        if not os.path.exists(self.config.annotate_csv):
            raise FileNotFoundError(f"Annotation file NOT FOUND: {self.config.annotate_csv}")
        self._check_csv_columns(
            self.config.annotate_csv,
            ANNOTATE_REQUIRED_COLS,
            "--annotate_file"
        )
        self.logger.info(f"Annotation file validated: {self.config.annotate_csv}")

        # bulkCNV.csv is optional — only validate if provided
        if self.config.bulk_csv:
            if not os.path.exists(self.config.bulk_csv):
                raise FileNotFoundError(f"Bulk CNV file NOT FOUND: {self.config.bulk_csv}")
            self._check_csv_columns(
                self.config.bulk_csv,
                BULK_CNV_REQUIRED_COLS,
                "--bulkCNV_file"
            )
            self.logger.info(f"Bulk CNV file validated: {self.config.bulk_csv}")
        else:
            self.logger.info("No bulk CNV file provided; proceeding without bulkCNV info.")

        wtcnr_dir = os.path.join(self.config.output_dir, "wtcnr")
        if not os.path.exists(wtcnr_dir):
            self.logger.warning(f"Weighted CNR directory {wtcnr_dir} not found. Step 6 may fail.")

    def _check_csv_columns(self, filepath, required_cols, arg_name):
        """Read only the header row and verify required column names are present."""
        actual_cols = set(pd.read_csv(filepath, nrows=0).columns)
        missing = required_cols - actual_cols
        if missing:
            raise ValueError(
                f"{arg_name} is missing required column(s): {missing}. "
                f"Expected columns: {required_cols}, got: {actual_cols}"
            )
