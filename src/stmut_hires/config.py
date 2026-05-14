import logging
from importlib import resources

logging.basicConfig(level=logging.INFO)


class BaseConfig:
    """Base configuration shared across all pipeline steps."""

    def __init__(self, output_dir, dry_run=False, cutoff=1000, window=100, num_processes=None, cores=4):
        self.dry_run = dry_run
        self.output_dir = output_dir
        self.cutoff = cutoff
        self.window = window
        self.num_processes = num_processes
        self.cores = cores
        self.logger = logging.getLogger(self.__class__.__name__)

        data_path = resources.files("stmut_hires.data")
        self.gene_info = str(data_path.joinpath("ensembl-gene-info.hg38.tsv"))
        self.corr_file = str(data_path.joinpath("tcga-skcm.cnv-expr-corr.tsv"))
        self.bed_file = str(data_path.joinpath("reference").joinpath("hg38_centromereSimple.bed"))


class InitialStepConfig(BaseConfig):
    """Configuration for step 1: reading cluster and expression inputs."""

    def __init__(
        self, 
        clusterf=None, 
        exp_h5=None, 
        bw_method=0.1,
        filter_cutoff=None,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.clusterf = clusterf
        self.exp_h5 = exp_h5
        self.bw_method = bw_method
        self.filter_cutoff = filter_cutoff
        
        self.expression_file = None
        self.ensembl_file = None

    def set_step1_outputs(self, expression_file_path, ensembl_file_path):
        """Update config with file paths produced by Step 1."""
        self.expression_file = expression_file_path
        self.ensembl_file = ensembl_file_path
        self.logger.info("Config updated with Step 1 outputs.")


class MergerConfig(InitialStepConfig):
    """Configuration for step 2: spatial barcode merging."""

    def __init__(self, spatial_file=None, **kwargs):
        super().__init__(**kwargs)
        self.spatial_file = spatial_file


class CNVAnalysisConfig(MergerConfig):
    """Configuration for step 6: CNV calling and visualization."""

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
