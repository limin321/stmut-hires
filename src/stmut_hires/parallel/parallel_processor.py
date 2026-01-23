import glob
import logging
import multiprocessing


from stmut_hires.parallel.spatial_loader import SpatialDataLoader
from stmut_hires.parallel.output_manager import OutputDirManager
from stmut_hires.parallel.cluster_merger import ClusterMerger

class ParallelClusterProcessor:
    """Responsible for parallel processing of clusters"""
    def __init__(self,config):
        self.config = config
        self.output_dir=config.output_dir
        self.spatial_file = config.spatial_file
        self.num_processes = getattr(config, 'num_processes', None) or multiprocessing.cpu_count() - 1
        self.cutoff = config.cutoff
        self.window = config.window
        self.logger = logging.getLogger(__name__)

    def parallel_cluster_merger(self):
        """Parallelly process clusters barcodes merging"""
        # Use the config's expression_file pattern which is already set in main.py
        clusters = glob.glob(self.config.expression_file)

        if not clusters:
            self.logger.error(f"No cluster parquet files found matching patter: {self.config.expression_file}")
            return # Exit if no files found
        
        # Load spatial_coordinates file
        spatial_loader = SpatialDataLoader(self.spatial_file)
        coords = spatial_loader.load_spatial_data()
        ensembl_path = self.config.ensembl_file

        # Create txt folder
        output_creator = OutputDirManager(self.output_dir)
        txt_output_dir, *others = output_creator.create_output_dir()

        with multiprocessing.Pool(processes=self.num_processes) as pool:
            args = [(cluster, coords, txt_output_dir, ensembl_path, self.cutoff, self.window) for cluster in clusters]
            pool.starmap(self._process_cluster_wrapper, args)

    @staticmethod
    def _process_cluster_wrapper(cluster_path, coords, txt_output_dir, ensembl_path, cutoff, window):
        """Wrapper method for multiprocessing"""
        merger = ClusterMerger(coords, txt_output_dir, ensembl_path, cutoff, window)  
        merger.process_cluster(cluster_path)      

