import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


#class CNVScoreQQPlot:
class AnalysisVisualizer:
    """
    Responsible for generating various analysis plot. 
    """
    def __init__(self, permuts_long, caseCNVscore_FDR, pmtimes, output_dir):
        self.permuts_long = permuts_long
        self.caseCNVscore_FDR = caseCNVscore_FDR
        self.pmtimes = pmtimes
        self.output_dir = output_dir
        self.output_dir_figures = os.path.join(self.output_dir, "figures")

    def qq_data(self):
        """Generate data for CNVScore QQ plot """
        pemts = self.permuts_long['permutScore']
        downS = []
        i = int(self.pmtimes/2)
        for n in range(pemts.shape[0]//self.pmtimes):
            a = pemts.iloc[i]
            downS.append(a)
            i = i + self.pmtimes

        data = pd.DataFrame({
            'CNVScore':self.caseCNVscore_FDR['CNVScore'].values, 
            'PermutScore':downS
        })
        return data
    
    def cnv_score_qq_plot(self, data):
        """Create QQ scatter plot of CNVscore and permutation values"""
        # Create the plot
        plt.figure(figsize=(8, 6))
        plt.scatter(data['PermutScore'], data['CNVScore'], color='lightgrey', edgecolor='black')

        # Add the y=x reference line
        plt.axline(xy1=(0, 0), slope=1, color='black', linestyle='--') # linestyle='--' is lty=2

        # Add labels and title (optional)
        plt.xlabel("PermutScore")
        plt.ylabel("CNVscore")
        plt.grid(True)
        plt.title("QQ Plot of CNV Scores vs. Permutations")
        plt.margins(0.08)

        output_path = os.path.join(self.output_dir_figures, "CNVs_RankedBySimilarityToDNA_QQplot.pdf")
        plt.savefig(output_path, format='pdf', bbox_inches='tight')
        plt.show()

    def plot_barcode_histogram(self, merged_df):
        """Generate a histogram of log-transformed merged barcode counts. """
        merged_df['log_barcode_count'] = np.log1p(merged_df['merged_barcodes_count'])
        plt.figure(figsize=(8, 5))
        plt.hist(
            x=merged_df['log_barcode_count'], 
            bins=200, 
            color='black',          
            edgecolor='grey'        
        )
        plt.title("Size Distribution of Merged Barcode Groups")
        plt.xlabel("Count of Constituent Barcodes (log1p scale)")
        plt.ylabel("Number of Merged Groups")
        plt.grid(False)
        
        output_path = os.path.join(self.output_dir_figures,"barcodes_counts_histogram.pdf")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.show() # or return the figure object

    def plot_cnvscore_histogram(self):
        """Generate Case CNVScore Histogram."""
        scores = self.caseCNVscore_FDR['CNVScore']
        plt.figure(figsize=(10,6))
        plt.hist(scores, bins=50, color='skyblue', edgecolor='black')
        plt.title('Histogram of CNV Scores')
        plt.xlabel('CNV Scores')
        plt.ylabel('Number of CNVScores')
        plt.savefig(os.path.join(self.output_dir_figures, "CNVs_RankedBySimilarityToDNA_CNVscoreHistogram.pdf"))
        plt.show()






