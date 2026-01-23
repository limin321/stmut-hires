import pandas as pd
import numpy as np

class ClusterTotalReads:
    """Responsible for Cluster TotalReads Data"""
    @staticmethod
    def _tumor_normal_label(sorted_reads):
        """Prepare tumor and normal labels"""
        sorted_cluster_tumor = sorted_reads[['cluster', 'annotate']].copy()
        sorted_cluster_tumor['label'] = 'T'
        sorted_cluster_tumor.loc[sorted_cluster_tumor['annotate']=='normal', 'label'] = 'N'
        sorted_cluster_tumor = sorted_cluster_tumor.drop_duplicates(subset=['cluster'], keep='first')
        return sorted_cluster_tumor
    
    @staticmethod
    def reads_his_params(cdt_sortedby_cluster_totalreads):
        """Prepare parameters for cluster total-reads histogram"""
        sorted_reads = cdt_sortedby_cluster_totalreads[['barcode','cluster', 'annotate','TotalRDs']]
        sorted_reads = sorted_reads.copy()
        sorted_reads['TotalRDs'] = pd.to_numeric(sorted_reads['TotalRDs'], errors='coerce')

        reads_normalized = []
        labels_reads = []
        sequential_lens = []

        uniq_clusters = sorted_reads['cluster'].unique()
        for cluster in uniq_clusters:
            clst_dat = sorted_reads[sorted_reads['cluster']==cluster]
            max_reads = clst_dat['TotalRDs'].max()
            if max_reads > 0:
                rate = clst_dat['TotalRDs']/max_reads
            else:
                rate = clst_dat['TotalRDs']

            reads_normalized.extend(rate.values)
            sequential_lens.append(len(rate))

            t1 = np.where(uniq_clusters == cluster)[0][0] + 1 # python 0-based, so need to +1 to start 1
            if cluster == uniq_clusters[0]:
                label_loc = 0.5 *len(rate)
                labels_reads.append(label_loc)
            else:
                label_loc = 0.5 *(len(rate)) + sum(sequential_lens[:t1-1])
                labels_reads.append(label_loc)

        sorted_cluster_tumor = ClusterTotalReads._tumor_normal_label(sorted_reads)
        return sorted_cluster_tumor,reads_normalized,labels_reads,uniq_clusters
