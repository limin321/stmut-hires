import os
import multiprocessing as mp
import logging 
from importlib import resources

# Assuming these imports are correct
from .cnr_reader import CnrReader
from .centromere_reader import CentromereReader
from .weighted_median import ArmWeightedMedian
from .cnr_writer import WtMedianWriter

class WorkflowOrchestrator:
    """Manage Overall Workflow,file discovery, parallel execution using existing class components"""
    def __init__(self, config):
        self.config = config
        self.output_dir = config.output_dir
        self.cores = config.cores
        self.cnr_input_dir = os.path.join(self.output_dir, "cnr")
        self.cnr_files = self._find_cnr_files()
        self.logger = logging.getLogger(__name__)

        # # Determine project root
        # module_dir = os.path.dirname(os.path.abspath(__file__))  # adjust if using Jupyter
        # project_root = os.path.abspath(os.path.join(module_dir, "../../.."))
        # self.bed_file = config.bed_file or os.path.join(project_root, "data/reference/hg38_centromereSimple.bed")
        self.bed_file = config.bed_file

        # check if bed_file exist
        if not os.path.exists(self.bed_file):
            error_msg = f"Required BED file not found: {self.bed_file}"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)

    def _find_cnr_files(self):
        """Discovers all CNR files in the input directory."""
        if not os.path.exists(self.cnr_input_dir):
            self.logger.warning(f"Input directory not found: {self.cnr_input_dir}")
            return []

        full_paths = [os.path.join(self.cnr_input_dir, f) 
                      for f in os.listdir(self.cnr_input_dir) 
                      if f.endswith('.cnr')]
        return full_paths

    @staticmethod
    def _single_cnr_weighted_median(cnr_file,bed_file,output_dir):
        """Calculate single cnr weighted median
        This function must be static/standalone for multiprocessing to work well.
        """
        try: 
            # 1. Read CNR file
            cnrloader = CnrReader(cnr_file)
            cnr_name, cnr = cnrloader.read_cnr()

            # 2. Read Centromere data (ideally once outside the loop if possible, 
            #    but kept here to keep the helper self-contained for multiprocessing simplicity)
            ctmereloader = CentromereReader(bed_file)
            centm_sorted = ctmereloader.read_centromere() 
            # 3. Calculate Arm median
            armsmedian = ArmWeightedMedian(cnr,centm_sorted)
            df_arms, df_wmedian = armsmedian.chr_arm_weighted_median()
            # 4. Write output
            savedata = WtMedianWriter(output_dir,df_arms,df_wmedian,cnr_name)
            savedata.save_wt_cnr()
            return f"[SUCCESS] Processed {cnr_name}"

        except Exception as e:
            return f"[ERROR] Failed to process {cnr_file}: {e}"


    def run_parallel_workflow(self):
        """Executes the single-file processing function in parallel."""
        if not self.cnr_files:
            self.logger.warning(f"No .cnr files found in the {self.cnr_input_dir} folder")
            return

        self.logger.info(f"Found {len(self.cnr_files)} cnr files. Starting parallel processing on {self.cores} cores")

        # Prepare argument tuples for starmap
        tasks = [(cnr_file, self.bed_file, self.output_dir) for cnr_file in self.cnr_files]
        with mp.Pool(processes=self.cores) as pool:
            results = pool.starmap(WorkflowOrchestrator._single_cnr_weighted_median, tasks)

   
