import logging
import os
import sys
from typing import Optional # version-safe, "str | None only for python 3.10+"

# new_path = "/stomics_data/liminData/Visium/stmut_python/stmut-hires/python/"
# sys.path.insert(0, new_path)

from call_cnv.barcodes_summary import ClusterSummary
from call_cnv.cluster_annotate import LoadAnnotate
from call_cnv.annotate_barcode import IntegrateAnnotate
from call_cnv.load_cdt import Loadcdt
from call_cnv.infer_cnv import InferCNV
from call_cnv.save_cnv import CNVWriter
from visualization.chromosome_meta import SortedChrom
from visualization.cluster_totalreads import ClusterTotalReads
from visualization.cnv_sorted_bycluster_plot import CNVSortedbyClusters
from visualization.unrooted_cnv_clustered_plot import UnrootedCNVHeatmap
from visualization.analysis_plots import AnalysisVisualizer
from call_cnv.cnvscore_orchestrators import CNVSortedByBulkCNV
from call_cnv.barcode_quintile import AssignBarcodeQuintile


class CNVCallPlotWorkflow:
    """Responsible for calling CNV and plot sorted by cluster
    CNVCallPlotWorkflow is responsible for performing CNV (Copy Number Variation) analysis 
    and generating related plots from cdt/grpWt.cdt.

    The workflow consists of several key steps:

    1. CNV Calling (`call_cnv`):
       - Summarizes barcodes by cluster.
       - Loads and processes the CDT expression matrix.
       - Normalizes CNV data for non-tumor references.
       - Sorts CNV by cluster and total reads.
       - Saves the merged barcode summary and sorted CNV data.

    2. CNV Plotting (`plot_cnv`):
       - Prepares chromosome metadata.
       - Generates cluster-sorted CNV heatmaps.
       - Generates unrooted CNV clustered heatmaps.

    3. Bulk CNV Analysis (optional, `plot_cnv_bulk_sorted`):
       - If bulk CNV data is provided, processes CNV scores against permutation results.
       - Generates QQ plots, histograms, and quintile assignments for barcodes.

    4. Full Workflow (`cnv_workflow`):
       - Runs CNV calling and plotting sequentially.
       - Executes optional bulk CNV workflow if data is provided.

    Attributes:
        graph_based_csv (str): Path to the CSV containing graph-based barcode clustering info.
        annotate_csv (str): Path to the annotation CSV for barcodes.
        output_dir (str): Directory for saving outputs (tables and plots).
        bulk_csv (str, optional): Path to optional bulk CNV file.
        pmtimes (int, optional): Number of permutations for bulk CNV analysis.
        merged_df (pd.DataFrame): DataFrame storing merged barcode summaries.
        cdt_meta: Metadata from CDT file.
        cdt_sortedby_cluster_totalreads: Sorted CNV data matrix.
        clusterSortedcnv: CNV data organized for plotting.
    """
    def __init__(self, graph_based_csv, annotate_csv, output_dir):
        """ 
        Initializes the CNVCallPlotWorkflow instance with input files and output directory.
        Optional bulk parameters can be set later if bulk data is available.
        """
        self.logger = logging.getLogger(__name__)
        self.graph_based_csv = graph_based_csv
        self.annotate_csv = annotate_csv
        self.output_dir = output_dir

        if not os.path.exists(self.graph_based_csv):
            raise FileNotFoundError(f"Missing graph_based_csv: {self.graph_based_csv}")

        if not os.path.exists(self.annotate_csv):
            raise FileNotFoundError(f"Missing annotate_csv: {self.annotate_csv}")

        self.dirpath = os.path.join(self.output_dir,"cluster_summary")
        self.cdt_file = os.path.join(self.output_dir, "cdt/grpWt.cdt")

        # results storage
        self.merged_df = None
        self.cdt_meta = None
        self.cdt_sortedby_cluster_totalreads = None
        self.clusterSortedcnv =None
    
        # Ensure output directories exist
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "tables"), exist_ok=True)

    def call_cnv(self):
        """Call CNV and save merged summary."""
        self.logger.info("Starting CNV calling...")
        clustersum = ClusterSummary(self.dirpath,self.graph_based_csv)
        merged_df = clustersum.merge_cluster_summary()
        merged_barcode_summary = clustersum.cluster_summary()
        annotate = LoadAnnotate.load_annotate(self.annotate_csv)       
        annotated_barcode = IntegrateAnnotate(annotate,merged_barcode_summary).annotate_barcode()

        # Load, normalize cdt, call CNV.
        cdt_matrix, cdt_meta = Loadcdt.cdt_loader(self.cdt_file)
        cdt = InferCNV(cdt_matrix,annotated_barcode)
        nontumor_normalized_cdt = cdt.cdt_processor() # same as c6 in R
        # sort cnv first by cluster then by total reads
        cdt_sortedby_cluster_totalreads = InferCNV.sort_cluster_totalreads(annotated_barcode,nontumor_normalized_cdt)
        cluster_sorted_cnv = CNVWriter.save_cluster_sorted_cnv(cdt_sortedby_cluster_totalreads, cdt_meta, self.output_dir)
        
        # Store these results as instance variables if they are needed by other methods(like plot_cnv)
        self.merged_df = merged_df
        self.cdt_meta = cdt_meta
        self.cdt_sortedby_cluster_totalreads = cdt_sortedby_cluster_totalreads
        self.clusterSortedcnv = cluster_sorted_cnv
        CNVWriter.save_merged_df(merged_df, self.output_dir)
        return merged_df,cdt_meta,cdt_sortedby_cluster_totalreads,cluster_sorted_cnv

    def plot_cnv(
        self, 
        ncluster=8,
        distance_metric: Optional[str] = None,
        linkage_method: Optional[str] = None,
    ):
        """ Prepare chromosome and cluster, totalreads for plotting."""
        self.logger.info(f"Ploting CNV heatmap sorted by cluster, total_reads, and unrooted hierarchical.")
        chr_meta = SortedChrom.gene_count_per_chr(self.cdt_meta)
        (sorted_cluster_tumor,reads_normalized,labels_reads,uniq_clusters) = ClusterTotalReads.reads_his_params(self.cdt_sortedby_cluster_totalreads)
        
        # --1. plot cluster and total reads sorted cnv heatmap.--
        CNVSortedbyClusters.cluster_sorted_cnv_heatmap(
            self.clusterSortedcnv,
            reads_normalized,
            labels_reads,
            sorted_cluster_tumor,
            uniq_clusters,
            chr_meta,
            self.output_dir
        )

        # --2. Plot unrooted_clustered_CNV_heatmap. --
        heatmap_gen = UnrootedCNVHeatmap(clusterSortedcnv=self.clusterSortedcnv, chr_meta=chr_meta, output_dir=self.output_dir)
        heatmap_gen.set_clustering_parameters(
            ncluster=ncluster,
            distance_metric=distance_metric,
            linkage_method=linkage_method,
        )
        # Generate the heatmap
        heatmap_gen.generate_heatmap(output_file_name=f"unrooted_CNVs_clustered_heatmap_class_{ncluster}clusters.pdf")
        self.logger.info("CNV plotting completed.")

    def plot_cnv_bulk_sorted(self,permuts_long,caseCNVscore_FDR, pmtimes):
        """Generate bulk CNV plots including QQ plot, histogram, and quintiles. """
        self.logger.info("Plotting bulk_sorted CNV results...")
        visualizer = AnalysisVisualizer(permuts_long, caseCNVscore_FDR, pmtimes, self.output_dir)
        data = visualizer.qq_data()
        visualizer.cnv_score_qq_plot(data)
        visualizer.plot_cnvscore_histogram()  
        visualizer.plot_barcode_histogram(self.merged_df)

        # Generate quintile data
        all_barcodes_quintiles=AssignBarcodeQuintile.barcode_quintils(self.merged_df,caseCNVscore_FDR)
        CNVWriter.save_qt(all_barcodes_quintiles, self.output_dir)
        self.logger.info("CNV bulk-sorted plotting completed.")

    def cnv_workflow(
        self, 
        ncluster=13, 
        bulk_csv=None, 
        pmtimes=None,
        distance_metric=None,
        linkage_method=None,
        ):
        """
        Run the full CNV workflow.

        Args:
            ncluster (int): Number of clusters for unrooted heatmap.
            bulk_csv (str, optional): Path to bulk CNV CSV for optional analysis.
            pmtimes (int, optional): Number of permutations for bulk CNV.
        """
        self.call_cnv()
        self.plot_cnv(
            ncluster=ncluster,
            distance_metric=distance_metric,
            linkage_method=linkage_method,
        )

        if bulk_csv is not None:
            if pmtimes is None:
                raise ValueError("pmtimes must be set when bulk_csv is provided")

            cnvSortedBulk = CNVSortedByBulkCNV(
            bulk_csv, 
            pmtimes, 
            self.output_dir, 
            self.clusterSortedcnv
            )
            permuts_long, caseCNVscore_FDR = cnvSortedBulk.run_pipeline()
            self.plot_cnv_bulk_sorted(permuts_long, caseCNVscore_FDR, pmtimes)
        else:
            self.logger.info("No bulk data provided, skipping bulk CNV workflow.")



