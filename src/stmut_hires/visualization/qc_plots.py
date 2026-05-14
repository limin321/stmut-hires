import os
import matplotlib.pyplot as plt


class QCplot:
    """ 
    Plot Histomgram of gene counts per cell before and after filtering.
    """
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.output_dir_figures = os.path.join(self.output_dir,"figures")


    def qc_hist(self, df1, df2):

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

        # Before filtering
        axes[0].hist(
            df1["gene_counts"], 
            bins=200, 
            color='skyblue',
            edgecolor='black',
            alpha=0.7
        )

        axes[0].set_xlabel("Gene Counts")
        axes[0].set_ylabel('Count')
        axes[0].set_title("Before Filtering")

        # After filtering
        axes[1].hist(
            df2["gene_counts"], 
            bins=200, 
            color='skyblue',
            edgecolor='black',
            alpha=0.7
        )

        axes[1].set_xlabel("Gene Counts")
        axes[1].set_title("After Filtering")

        plt.tight_layout()
        os.makedirs(self.output_dir_figures, exist_ok=True)
        output_path = os.path.join(self.output_dir_figures, "gene_counts_before_after_filtering.png")
        plt.savefig(output_path, dpi=300,bbox_inches="tight")
        plt.show()