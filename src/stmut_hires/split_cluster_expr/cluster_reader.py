# Handle cluster read and process
import pandas as pd

class ClusterReader:
    """Read and process cluster information """
    @staticmethod
    def read_clusters(cluster_file):
        cluster_data = pd.read_csv(cluster_file)
        cluster_data['Graph-based'] = cluster_data['Graph-based'].str.replace(" ","")

        return {
            cluster: cluster_data.loc[cluster_data['Graph-based']==cluster, "Barcode"].tolist()
            for cluster in cluster_data["Graph-based"].unique()
        }
