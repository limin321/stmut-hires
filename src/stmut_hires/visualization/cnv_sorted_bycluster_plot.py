import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
import seaborn as sns

class CNVSortedbyClusters:
    """Responsible for plotting cluster sorted CNV heatmap"""
    
    @staticmethod
    def _cnv_params(clusterSortedcnv):
        # 1.1 Prepare cnv_matrix for heatmap
        cnv_matrix = clusterSortedcnv.drop(['CLID', 'NAME'], axis=1).copy()
        cnv_data = cnv_matrix.values # Get the underlying numpy array for efficient plotting
        cnv_data = cnv_data[::-1]

        # 1.2 Define the colors: Set vmin and vmax for a data-driven color scale
        #vmin_val = np.min(cnv_matrix.values)
        #vmax_val = np.max(cnv_matrix.values)
        norm = Normalize(vmin=-1, vmax=0.8) # If you are using Matplotlib's pcolormesh or imshow directly, use this:

        myCol = sns.diverging_palette(
            h_neg=240,  # Blue hues
            h_pos=20,   # Red/Orange hues
            l=80,       # Lightness: 70 is light/muted, 50 is default
            s=80,       # Saturation
            center='light', # Ensure white center
            as_cmap=True
        )
        return cnv_data, norm, myCol
        
    @staticmethod
    def cluster_sorted_cnv_heatmap(
        clusterSortedcnv,
        reads_normalized,
        labels_reads,
        sorted_cluster_tumor,
        uniq_clusters,
        chr_meta,
        output_dir
    ):
        """Plot CNV sorted by cluster and total reads"""
        # Create a new blank figure object - the overall container for all plots
        fig = plt.figure(figsize=(10, 10))

        # Creates the grid layout manager. 2,2-specifies a 2-row by 2-col grid; 
        gs = GridSpec(2, 2, figure=fig, width_ratios=[0.2, 2], height_ratios=[0.1, 2],
                    hspace=0.001, wspace=0.001) # hspace; wspace -- space between plots horizontally and vertically
        # Define axes:
        ax_empty = fig.add_subplot(gs[0, 0])      # Top Left (Empty)
        ax_cluster = fig.add_subplot(gs[0, 1])    # Top Right (Cluster Plot)
        ax_chr_bar = fig.add_subplot(gs[1, 0])    # Bottom Left (Chromosome Barplot)
        ax_cnv_image = fig.add_subplot(gs[1, 1])  # Bottom Right (Main CNV Heatmap)
        # Hide the empty axis
        ax_empty.axis('off')

        # --- 1. Plotting CNV Image (Bottom Right Panel) ---
        #cnv_data, norm, myCol = ClusterSortedPlot._cnv_params(clusterSortedcnv)
        cnv_data, norm, myCol = CNVSortedbyClusters._cnv_params(clusterSortedcnv)

        #  visualize cnv
        ax_cnv_image.imshow(cnv_data, cmap=myCol, norm=norm, aspect='auto', origin='lower', interpolation='nearest')
        ax_cnv_image.set_xlabel("Barcodes", fontsize=10)
        ax_cnv_image.set_xticks([])
        ax_cnv_image.set_yticks([])
        ax_cnv_image.set_title("") # Mimics the specific '20' label position

        # Remove the black box border around the heatmap 
        ax_cnv_image.spines['top'].set_visible(False)
        ax_cnv_image.spines['right'].set_visible(False)
        ax_cnv_image.spines['bottom'].set_visible(False)
        ax_cnv_image.spines['left'].set_visible(False)

        # --- 2. Plotting Cluster Data (Top Right Panel) ---
        colors1 = ["#929292"] * len(reads_normalized)
        ax_cluster.fill_between(range(1, len(reads_normalized) + 1), reads_normalized, 0, color=colors1, linewidth=0)

        # Remove internal padding:
        ax_cluster.set_xlim(1, len(reads_normalized) + 1) # Set x-limits exactly to data bounds
        ax_cluster.margins(x=0, y=0)                      # Set subplot margins to zero
        ax_cluster.set_title("Cluster", loc='center', pad=-10, fontsize=10)
        ax_cluster.axis('off')

        # Add Cluster/T/N labels as text
        for idx, loc in enumerate(labels_reads):
            y_pos_tn = 0.8 if idx % 2 == 0 else 0.6 # Stagger Y position for every second label
            y_pos_cluster = 0.5 if idx %2 == 0 else 0.3
            ax_cluster.text(loc, y_pos_tn, sorted_cluster_tumor['label'].iloc[idx], ha='center', va='top', fontsize=6, color='Black')
            ax_cluster.text(loc, y_pos_cluster, uniq_clusters[idx].replace('cluster',''), ha='center', va='center', fontsize=6, color='Black')

        # --- 3. Plotting Chromosome Data (Bottom Left Panel) ---
        COLOR_LIGHT_GREY = '#D3D3D3'
        COLOR_DARK_GREY = '#A9A9A9'
        TWO_COLORS = [COLOR_LIGHT_GREY, COLOR_DARK_GREY]
        y1 = chr_meta['count'].values # 'chr_meta' is pandas DataFrame and 'y1' is the list/array of counts
        chromosome_labels = chr_meta['chr'].values # e.g., [1, 2, 3, ...]

        # Calculate the cumulative sum of y1 values to determine the 'bottom' position of each bar/chr segment
        y_bottoms = np.cumsum(y1) - y1 
        total_height = np.sum(y1)
        ax_chr_bar.clear()

        for i, value in enumerate(y1):
            segment_color = TWO_COLORS[i % 2] 
            ax_chr_bar.bar(
                x=0.5,                  
                height=value,           
                width=1.0,              
                bottom=y_bottoms[i],    
                color=segment_color,    
                edgecolor='none',       
                align='center'
            ) # Plot the bar segment
            
            # --- Add the Chromosome Label inside the bar ---
            center_y = y_bottoms[i] + (value / 2)
            label_text = str(chromosome_labels[i]) # Get the label from your meta dataframe
            
            ax_chr_bar.text(
                x=0.5,             # X position: Center of the bar (0.5)
                y=center_y,        # Y position: Center of the current segment
                s=label_text,      # The label itself (e.g., '1', '2', '3')
                ha='center',       # Horizontal alignment: Centered
                va='center',       # Vertical alignment: Centered
                fontsize=8,        # Adjust font size as needed
                color='black'      # Set text color to black for visibility
            )

        # Set the Y-axis limit to the total height of all segments
        ax_chr_bar.set_ylim(0, total_height)
        ax_chr_bar.set_xlim(0, 1) # Set consistent X-axis limits to center the bar
        ax_chr_bar.set_ylabel("Chromosome", fontsize=10, loc='center')
        # Manually remove the spines (borders) and ticks 
        ax_chr_bar.spines['top'].set_visible(False)
        ax_chr_bar.spines['right'].set_visible(False)
        ax_chr_bar.spines['bottom'].set_visible(False)
        ax_chr_bar.spines['left'].set_visible(False)
        ax_chr_bar.set_xticks([])
        ax_chr_bar.set_yticks([])
        ax_cluster.set_xticks([])
        ax_cluster.set_yticks([])

        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        figpath = f"{output_dir}/figures/CNVs_OrganizedByGEcluster_UMIcount.pdf"
        plt.savefig(figpath)
        plt.show()
