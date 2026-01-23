import pandas as pd
import numpy as np
import logging
import scipy.cluster.hierarchy as sch 
from scipy.spatial.distance import pdist
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

class UnrootedCNVHeatmap:
    """Responsible for generating a complex clustered heatmap for CNV data"""
    def __init__(self, clusterSortedcnv: pd.DataFrame, chr_meta: pd.DataFrame, output_dir: str):
        """ 
        Initialize the generator with pre-loaded Pandas DataFrame
        """
        self.logger = logging.getLogger(__name__)
        self.output_dir = output_dir
        self.labelCdt = clusterSortedcnv.iloc[:, 0:2]
        self.cnv_matrix = clusterSortedcnv.drop(['CLID', 'NAME'], axis=1)
        self.cnv_matrix_T = self.cnv_matrix.transpose()
        self.cnv_matrix_T = self.cnv_matrix_T.iloc[:500,] # for quick test
        self.chr_meta = chr_meta[::-1]
        self.ncluster = 6
        self.method1 = 'canberra'
        self.method2 = 'complete'
        self.g = None
        self.logger.info("Data initialized successfully from DataFrame.")

    def set_clustering_parameters(self, ncluster=None, distance_metric=None, linkage_method=None):
        """Dynamically set the parameters used for hierarchical clustering. """
        if ncluster is not None:
            self.ncluster = ncluster
        if distance_metric is not None:
            self.method1 = distance_metric
        if linkage_method is not None:
            self.method2 = linkage_method
        self.logger.info(f"Clustering parameters updated: Clusters={self.ncluster}, Metric='{self.method1}', Linkage='{self.method2}'")
    
    def generate_heatmap(self, output_file_name="Unrooted_CNVs_clustered_heatmap.pdf"):
        """Performs clustering and generates the final heatmap visualization. """
        dist_matrix = pdist(self.cnv_matrix_T, metric=self.method1)
        cluster_row = sch.linkage(dist_matrix, method=self.method2)
        gr_row = sch.fcluster(cluster_row, t=self.ncluster, criterion='maxclust')
        col1 = sns.color_palette('Set1', self.ncluster)
        row_colors_list = pd.Series(gr_row).map(lambda x: col1[x-1]).tolist()

        # Define custom colormap
        colors = ["#1984c5", "#22a7f0", "#a7d5ed", "white", "white", "#e1a692", "#de6e56", "#c23728"]
        my_cmap = LinearSegmentedColormap.from_list('custom_cnv', colors, N=101)

        # Generate Heatmap
        self.g = sns.clustermap(
            self.cnv_matrix_T,
            method=self.method2,
            metric=self.method1,
            cmap=my_cmap,
            col_cluster=False,
            row_colors=row_colors_list,
            cbar_kws={'label': 'CNV Level'},
            xticklabels=False,
            yticklabels=False,
            figsize=(12, 11),
            cbar_pos=[0.05, 0.85, 0.008, 0.07] 
        )
        self.g.ax_cbar.tick_params(labelsize=6)
        self.g.ax_cbar.yaxis.label.set_size(8)

        # Add the chromosome Bar Plot
        self._add_chromosome_bar_plot(self.chr_meta)

        # Final Adjustment and Saving
        self._finalize_plot_appearance(col1)
        self.g.savefig(f"{self.output_dir}/figures/{output_file_name}", bbox_inches='tight')
        plt.show()

    def _add_chromosome_bar_plot(self, chr_meta):
        """Add the chr barplot above the heatmap """
        COLOR_LIGHT_GREY = '#D3D3D3'
        COLOR_DARK_GREY = '#A9A9A9'
        TWO_COLORS = [COLOR_LIGHT_GREY, COLOR_DARK_GREY]
        y1 = chr_meta['count'].values
        chromosome_labels = chr_meta['chr'].values
        x_bottoms = np.cumsum(y1) - y1
        total_width = np.sum(y1)
        ax_heatmap = self.g.ax_heatmap
        ax_chr_bar_pos = ax_heatmap.get_position()
        new_bottom = ax_chr_bar_pos.y1 + 0.001
        height = 0.02
        ax_chr_bar = self.g.fig.add_axes([ax_chr_bar_pos.x0, new_bottom, ax_chr_bar_pos.width, height])
        for i, value in enumerate(y1):
            segment_color = TWO_COLORS[i % 2]
            ax_chr_bar.bar(
                x=x_bottoms[i], height=1.0, width=value, bottom=0.0, 
                color=segment_color, edgecolor='none', align='edge'
            )
            center_x = x_bottoms[i] + (value / 2)
            ax_chr_bar.text(x=center_x, y=0.5, s=str(chromosome_labels[i]),
                            ha='center', va='center', fontsize=8, color='black')
        ax_chr_bar.set_xlim(0, total_width)
        ax_chr_bar.set_ylim(0, 1)
        ax_chr_bar.set_xticks([])
        ax_chr_bar.set_yticks([])
        ax_chr_bar.spines['top'].set_visible(False)
        ax_chr_bar.spines['bottom'].set_visible(False)
        ax_chr_bar.spines['right'].set_visible(False)
        ax_chr_bar.spines['left'].set_visible(False)
        ax_chr_bar.set_title("Chromosomes", fontsize=10, pad=5)
        self.logger.info("Chromosome bar plot added.")

    def _finalize_plot_appearance(self, col1):
        """Helper method for titles, labels, and the manual legend."""
        self.g.fig.suptitle("Copy Number Variations", fontsize=16, y=0.9)
        self.g.ax_heatmap.set_title("")
        self.g.ax_heatmap.set_ylabel("Barcodes")
        self.g.ax_heatmap.set_xlabel("Genes")
        for label in range(1, self.ncluster + 1):
            self.g.ax_heatmap.bar(0, 0, color=col1[label-1], label=f'Cluster {label}', linewidth=0)
        self.g.ax_heatmap.legend(loc='upper left', bbox_to_anchor=(-0.1, 1.2), frameon=False, fontsize=8)
        self.logger.info("Heatmap appearance finalized.")

