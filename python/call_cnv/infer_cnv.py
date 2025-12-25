import pandas as pd
class InferCNV:
    """Rresponsible for inferring CNV from raw cdt file"""
    def __init__(self, cdt, annotate):
        self.cdt = cdt
        self.annotate = annotate

    def cdt_processor(self):
        """Normalize cdt file"""
        cdt1 = self.cdt.iloc[2:, 4:]
        cdt1 = cdt1.apply(pd.to_numeric, errors='coerce')

        # cap CNV values to [-1,1]
        cdt1[cdt1>1] = 1
        cdt1[cdt1 < -1] = -1
        c1 = cdt1.transpose()
        c2 = c1.reset_index(names='barcode')
        c3 = pd.merge(self.annotate[['barcode','cluster','annotate']],c2, on='barcode')
        c3 = c3.set_index('barcode')
        c4 = c3.iloc[:, 2:]
        c4 = c4.apply(pd.to_numeric, errors='coerce')

        # Calculate the proportion of zeros in each row/barcode
        zero_proportion = (c4==0).sum(axis=1)/(c4.shape[1] - 2)

        # keep rows where proportion of zeros is less than 0.05
        c4 = c4[zero_proportion < 0.05]

        # extract non-tumor/normal barcodes and calculate median of each gene
        nontumor_barcode_list = c3[c3['annotate'] == 'normal'].index.values
        non_tumor = c4[c4.index.isin(nontumor_barcode_list)]
        gene_median = non_tumor.median(axis=0)

        # subtract the median from each spot (row-wise operation)
        c4_nontumor_normalize = c4.subtract(gene_median, axis=1)
        #c4_nontumor_normalize = c4_nontumor_normalize.transpose()
        return c4_nontumor_normalize

    @staticmethod
    def sort_cluster_totalreads(annotated_barcode,c4_nontumor_normalize):
        cdt4sort=pd.merge(annotated_barcode,c4_nontumor_normalize, right_index=True, left_on='barcode', how='right')
        cdt4sort['TotalRDs'] = pd.to_numeric(cdt4sort['TotalRDs'], errors='coerce')
        cdt_sorted = cdt4sort.sort_values(
            by = ['cluster', 'TotalRDs'],
            ascending = [True, False],
            key = lambda x: x.str.lstrip('cluster').astype(int) if x.name == 'cluster' else x
        )
        return cdt_sorted
