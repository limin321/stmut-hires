import pandas as pd
import os


class ClusterSummary:
    """Combine the grouping info of individual clusters into one"""
    def __init__(self, dirpath, graph_based_csv):
        self.dirpath = dirpath
        self.graph_based_csv = graph_based_csv
        self.merged_df = None

    def merge_cluster_summary(self):
        files = os.listdir(self.dirpath)
        if not files:
            print("No Cluster*_barcode_grouping_info.csv found in directory. Move these files from txt dir to cluster_summary")
            return None

        f1 = f"{self.dirpath}/{files[0]}"
        merged_df = pd.read_csv(f1)
        for f in files[1:]:
            fpath = f"{self.dirpath}/{f}"
            df = pd.read_csv(fpath)
            merged_df = pd.concat([merged_df, df], ignore_index=True)
        
        # Count number of barcode in each merged_group
        merged_df['merged_barcodes_count'] = merged_df['merged_barcodes'].str.split(',').str.len().fillna(0).astype(int)
        self.merged_df = merged_df
        return merged_df

    def cluster_summary(self):
        """Generate a summary of merged barcodes with spatial data."""
        merged_df_sub = self.merged_df[['main_barcode','total_non_zero','array_row','array_col']]

        #read-in spatial location
        dat1 = pd.read_csv(self.graph_based_csv)
        dat = pd.merge(dat1, merged_df_sub, left_on='Barcode', right_on='main_barcode')
        dat = dat.drop('main_barcode', axis=1)
        merged_barcode_summary = dat[['Barcode','Graph-based','total_non_zero','array_row','array_col']]
        merged_barcode_summary.rename(columns={'Barcode':'barcode',
                            'Graph-based':'cluster',
                            'total_non_zero':'TotalRDs'}, inplace = True)
        merged_barcode_summary['cluster']=merged_barcode_summary['cluster'].str.replace(' ','').str.replace('Cluster', 'cluster')
        return merged_barcode_summary





