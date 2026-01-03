import os
import glob
import logging
import pandas as pd
from abc import ABC, abstractmethod
from importlib import resources

""" 
Use `super` to call the parent class
"""

# Set up basic logging
logging.basicConfig(level=logging.INFO)

class BaseConfig(ABC):
    """Base class for all configuration."""
    def __init__(self, output_dir, dry_run=False, cutoff=1000, window=100, num_processes=None, cores=4):
        self.dry_run = dry_run
        self._base_output_dir = output_dir #Store the base path
        self.output_dir = output_dir
        self.cutoff = cutoff
        self.window = window
        self.num_processes = num_processes
        self.cores = cores # for weighted-median parallel
        self.logger = logging.getLogger(self.__class__.__name__)

        # Get the path to the data resource
        data_pkg = "stmut_hires.data"
        # Use .files() to get a traversable path (Python 3.9+)
        data_path = resources.files(data_pkg)
        self.gene_info = str(data_path.joinpath("ensembl-gene-info.hg38.tsv"))
        self.corr_file = str(data_path.joinpath("tcga-skcm.cnv-expr-corr.tsv"))
        default_bed = str(data_path.joinpath("reference").joinpath("hg38_centromereSimple.bed"))
        self.bed_file = default_bed

    
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
    def __init__(self, clusterf=None, exp_h5=None, **kwargs):
        super().__init__(**kwargs)
        self.clusterf = clusterf
        self.exp_h5 = exp_h5
        self.expression_file = None # This will be set by step1 output
        self.ensembl_file = None    # This will be set by step1 output

    def validate_inputs(self):
        """Validate inputs for the first step."""
        self.logger.info("Validating step 1 input files...")
        if self.clusterf is not None:
            if not os.path.exists(self.clusterf):
                raise FileNotFoundError(f"Cluster file NOT FOUND: {self.clusterf}")
        
        if self.exp_h5 is not None:
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
        self.spatial_file = spatial_file
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
        if self.spatial_file is not None:
            if not os.path.exists(self.spatial_file):
                raise FileNotFoundError(f"Spatial input file not found: {self.spatial_file}")
        
        # We need the expression and ensembl files generated from the first step to exist
        # We can check if clusterf was provided as a proxy for "is this a full run?"
        if self.clusterf is not None:
            if not self.expression_file:
                raise FileNotFoundError(f"Expression file pattern is None: {self.expression_file}")
        
            files_found = glob.glob(self.expression_file)
            if not files_found:
                raise FileNotFoundError(f"No expression files found matching pattern: {self.expression_file}")
             
            if not self.ensembl_file or not os.path.exists(self.ensembl_file):
                raise FileNotFoundError(f"Ensembl file missing or not found after Step 1: {self.ensembl_file}")
            if self.spatial_file and not os.path.exists(self.spatial_file):
                raise FileNotFoundError(f"Spatial input file not found: {self.spatial_file}")
            
            # Validate general type constrains
            if not isinstance(self.cutoff, int):
                raise TypeError(f"cutoff must be an integer, but received: {type(self.cutoff)}")
            if not isinstance(self.window, int):
                raise TypeError(f"window must be an integer, but received: {type(self.window)}")

            self.logger.info("Step 2 merging barcodes inputs are valid.")

class CNVAnalysisConfig(MergerConfig):
    """Configuration for call_cnv and visualization modules."""
    def __init__(self, annotate_csv=None, bulk_csv=None, 
                 pmtimes=5, ncluster=6, distance_metric='euclidean', 
                 linkage_method='ward', **kwargs):
        super().__init__(**kwargs)
        self.annotate_csv = annotate_csv
        self.bulk_csv = bulk_csv
        self.pmtimes = pmtimes
        self.ncluster = ncluster
        self.distance_metric = distance_metric
        self.linkage_method = linkage_method


    def validate_inputs(self):

        # We only call super().validate_inputs() if we are doing a full 'run'
        # If clusterf is None, we are likely in 'call-cnv' mode
        if self.clusterf is not None:
            super().validate_inputs()
        else:
            self.ensure_output_dir()
            self.logger.info("Skipping Step 1-5 validation for call-cnv mode.")

        # Step 6 Validate: Annotation File
        if self.annotate_csv:
            if os.path.exists(self.annotate_csv):
                self.logger.info(f"Annotation file provided: {self.annotate_csv}")
            else:
                self.annotate_csv = None
        else:
            self.logger.info("No annotation file provided; using default cluster IDs.")

        # Validate Optional: Bulk CNV File
        if self.bulk_csv:
            if os.path.exists(self.bulk_csv):
                self.logger.info(f"Bulk CNV reference provided: {self.bulk_csv}")
            else:
                self.bulk_csv = None
        else:
            self.logger.info("No bulk CNV and annotation file provided; proceeding with standard workflow without bulkCNV info.")

        # Add check: Do Step 5 outputs exist?
        # Step 6 requires files in the output_dir created by previous runs
        wtcnr_dir = os.path.join(self.output_dir, "wtcnr") # adjust name to match your code
        if not os.path.exists(wtcnr_dir):
            self.logger.warning(f"Weighted CNR directory {wtcnr_dir} not found. Step 6 may fail.")