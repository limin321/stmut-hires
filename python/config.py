import os
import glob
import logging
import pandas as pd
from abc import ABC, abstractmethod

# Set up basic logging
logging.basicConfig(level=logging.INFO)

class BaseConfig(ABC):
    """Base class for all configuration."""
    def __init__(self, output_dir, dry_run=False, cutoff=1000, window=100):
        self.dry_run = dry_run
        self._base_output_dir = output_dir #Store the base path
        self.output_dir = output_dir
        self.cutoff = cutoff
        self.window = window
        self.num_processes = num_processes
        self.logger = logging.getLogger(self.__class__.__name__)
    
    @abstractmethod
    def validate_inputs(self):
        """Abstract method to be implemented by derived classes for specific validations."""
        pass

    def ensure_output_dir(self):
        """Helper to create the main output dir if not dry-run."""
        if not self.dry_run:
            # Check if the output_dir is actually a directory path, not just a file path
            if self.output_dir and not os.path.isdir(os.path.dirname(self.output_dir)) and os.path.dirname(self.output_dir) != '':
                pass
            os.makedirs(self.output_dir, exist_ok=True)
            self.logger.info(f"Created output dir: {self.output_dir}")
        else:
            self.logger.info(f"[DRY-RUN] Would create output dir: {self.output_dir}")


class InitialStepConfig(BaseConfig):
    """COnfiguration specific to the initial data processing step (step1)."""
    def __init__(self, clusterf, exp_h5, **kwargs):
        super().__init__(**kwargs)
        self.clusterf = clusterf
        self.exp_h5 = exp_h5
        self.expression_file = None # This will be set by step1 output
        self.ensembl_file = None    # This will be set by step1 output

    def validate_inputs(self):
        """Validate inputs for the first step."""
        self.logger.info("Validating step 1 input files...")
        if not os.path.exists(self.clusterf):
            raise FileNotFoundError(f"Cluster file NOT FOUND: {self.clusterf}")
        if not os.path.exists(self.exp_h5):
            raise FileNotFoundError(f"H5 file NOT FOUND: {self.exp_h5}")
        self.ensure_output_dir()
        self.logger.info("Step 1 inputs are valid.")

    def set_step1_outputs(self, expression_file_path, ensembl_file_path):
        """Method to upate config with paths generated during Step1."""
        self.expression_file = expression_file_path
        self.ensembl_file = ensembl_file_path
        self.logger.info(f"Config updated with Step 1 outputs.")

class MergerConfig(InitialStepConfig):
    """Configuration for the barcode merger step."""
    def __init__(self, spatial_file=None, **kwargs):
        super().__init__(**kwargs)
        self.spatial_file_path = spatial_file
        #self._spatial_df_loaded = None

    # @property
    # def spatial_df(self):
    #     #...(getter/setter methods as describing previously) ...
    #     if self._spatial_df_loaded is None:
    #         self.logger.warning("Accessing spatial_df before it is loaded!")
    #     return self._spatial_df_loaded
    
    # @spatial_df.setter
    # def spatial_df(self, df):
    #     if not isinstance(df, pd.DataFrame):
    #         raise TypeError(f"spatial_df must be a pandas DataFrame, but received: {type(df)}")
    #     self._spatial_df_loaded = df


    def validate_inputs(self):
        """Validate inputs for step2..."""
        # First, validate the inputs needed from the initial step
        super().validate_inputs()

        #Then, validate inputs specific to the merge step
        self.logger.info("Validating Step 2 specific inputs...")
        
        # We need the expression and ensembl files generated from the first step to exist
        if not self.expression_file:
             raise FileNotFoundError(f"Expression file pattern is None: {self.expression_file}")
        
        # --- MODIFIED: Use glob.glob() to check for existence of files matching the pattern ---
        files_found = glob.glob(self.expression_file)
        if not files_found:
             raise FileNotFoundError(f"No expression files found matching pattern: {self.expression_file}")
        # --- END MODIFIED ---
             
        if not self.ensembl_file or not os.path.exists(self.ensembl_file):
             raise FileNotFoundError(f"Ensembl file missing or not found after Step 1: {self.ensembl_file}")
        if self.spatial_file_path and not os.path.exists(self.spatial_file_path):
            raise FileNotFoundError(f"Spatial input file not found: {self.spatial_file_path}")
        
        # Validate general type constrains
        if not isinstance(self.cutoff, int):
            raise TypeError(f"cutoff must be an integer, but received: {type(self.cutoff)}")
        if not isinstance(self.window, int):
            raise TypeError(f"window must be an integer, but received: {type(self.window)}")

        self.logger.info("Step 2 merging barcodes inputs are valid.")



