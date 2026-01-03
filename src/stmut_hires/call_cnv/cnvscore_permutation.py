import pandas as pd
import numpy as np

class CNVScore:
    """Responsible for calculating CNV score, and permuation of cnv_scores"""
    @staticmethod
    def count_arm_genes(clusterSortedcnv):
        """Count the number of genes of each arm."""
        chromosome = [str(i) for i in range(1,25)]
        start = [122026459, 92188145, 90772458, 49712061, 46485900, 58553888, 58169653,
                    44033744, 43389635, 39686682, 51078348, 34769407, 16000000, 16000000,
                    17083673, 36311158, 22813679, 15460899, 24498980, 26436232, 10864560,
                    12954788, 58605579, 10316944]
        end = [124932724, 94090557, 93655574, 51743951, 50059807, 59829934, 61528020,
                45877265, 45518558, 41593521, 54425074, 37185252, 18051248, 18173523,
                19725254, 38265669, 26616164, 20861206, 27190874, 30038348, 12915808,
                15054318, 62412542, 10544039]

        centmere = pd.DataFrame({
            'chromosome': chromosome,
            'start': start,
            'end': end
        })

        # Extract chr, pos, gene from NAME column
        df1 = clusterSortedcnv[['NAME']].copy()
        df1[['chr', 'pos', 'gene']] = df1['NAME'].str.split(':', n=2, expand=True)

        # Split position into start and end
        df1[['start', 'end']] = df1['pos'].str.split('-', expand=True)
        df1['start'] = pd.to_numeric(df1['start'])
        df1['end'] = pd.to_numeric(df1['end'])
        df1 = df1[['NAME', 'chr', 'start', 'end']]
        df1['arms'] = ''

        # Iterate through chromosomes to assign p/q arms based on centromere location
        # Using .apply is more "pandas-like" than iterating row-by-row, but iterating
        # by chromosome groups is necessary to use the specific centromere start/end values.
        d1_new_list = []
        for ch in centmere['chromosome'].unique():
            d1 = df1[df1['chr'] == ch].copy()
            cent_start = centmere[centmere['chromosome']==ch]['start'].iloc[0]

            # Assign arms: if gene end < cent_start, it's 'p', otherwise 'q'
            d1['arms'] = np.where(d1['end'] < cent_start, f'{ch}p', f'{ch}q')
            d1_new_list.append(d1)

        # combine the results back into a single dataframe, maintaining original order
        d1_new = pd.concat(d1_new_list).sort_index()

        # Merge the new 'arms' column back into the original 'cdt' dataframe
        d1 = d1_new[['NAME', 'arms']]
        cdt = pd.merge(d1, clusterSortedcnv, on='NAME', how='inner')

        gene_counts = cdt.groupby('arms').size().reset_index(name='count')
        gene_counts['chr'] = gene_counts['arms'].str[:-1]
        gene_counts['chr'] = pd.to_numeric(gene_counts['chr'])
        gene_counts = gene_counts.sort_values(
            by = ['chr', 'arms'],
            ascending = [True, True]
        )

        first_row_indices = d1_new.groupby('arms').apply(lambda x: x.index.min())
        gene_counts = gene_counts.set_index('arms')
        gene_counts['gene_row'] = gene_counts.index.map(first_row_indices)
        gene_counts = gene_counts.reset_index()
        return gene_counts

    @staticmethod
    def _cal_cnv_score(df1_bulk, cdt, row_num_np):
        scores = []
        for i in range(cdt.shape[1]):
            df1_bulk['cell_value'] = cdt.iloc[row_num_np,i].values
            val = sum(df1_bulk['gainloss']*df1_bulk['wt']*df1_bulk['cell_value'])
            scores.append(val)
        return scores

    @staticmethod
    def sort_cnv_score(clusterSortedcnv,bulk,gene_counts,pmtimes):
        """Sorte cdt file based on cnv score."""
        name = clusterSortedcnv.iloc[:,:2].copy()
        cdt = clusterSortedcnv.iloc[:,2:].copy()

        df1 = gene_counts[gene_counts['arms'].isin(bulk['arms'])]
        df1_bulk = pd.merge(df1, bulk, on='arms', how='inner')

        # calculate weighted value of each arm
        genes = df1_bulk['count']
        max_gene_count = max(genes)
        row_num_np = np.array(df1_bulk['gene_row'])
        df1_bulk['wt'] = [n/max_gene_count for n in genes]
        cnv_scores = CNVScore._cal_cnv_score(df1_bulk, cdt, row_num_np)
        scores_indexed = pd.Series(cnv_scores, cdt.columns, name='CNVScore')

        # Sort cdt by CNVScore
        sorted_sample_columns = scores_indexed.sort_values(ascending=True).index
        cdt_sorted = cdt[sorted_sample_columns]

        # Case CNV score
        cdt_cnv_score_sorted = pd.concat([name, cdt_sorted], axis=1)
        caseCNVscore = scores_indexed.to_frame()
        caseCNVscore = caseCNVscore.reset_index(names=['barcode'])
        caseCNVscore.sort_values(by='CNVScore', ascending=False, inplace=True)
        caseCNVscore['Rank'] = caseCNVscore['CNVScore'].rank(method='dense', ascending=False).astype(int)

        # Permuted CNV scores
        permuted_scores_df = CNVScore.permut_scores(pmtimes,cdt_sorted,df1_bulk,row_num_np)
        return cdt_cnv_score_sorted, caseCNVscore, permuted_scores_df

    # Permutation
    @staticmethod
    def get_single_permut_cnv_score(cdt_sorted,df1_bulk,row_num_np):
        """Calculate permuted cdt cnv score."""
        arr = cdt_sorted.values # Get the underlying numpy array
        N_cols = arr.shape[1] # Get the number of columns (N_samples)

        # -- 1. Permut columns within each row independently --
        # For each row, generate a unique random permutation of column indices [0, 1, ..., N_cols-1]
        # np.random.rand returns an array of random floats (N_rows x N_cols)
        # np.argsort sorts these floats, giving the random order for each row
        row_permutations = np.argsort(np.random.rand(*arr.shape), axis=1)
        arr_permuted_cols = arr[np.arange(arr.shape[0])[:, None], row_permutations]

        # --2. Permute rows within each column independently --
        # This process is similar to Step 1, but applied vertically (axis=0)
        N_rows = arr.shape[0]

        # For each column, generate a unique random permutation of row indices [0, 1, ..., N_rows-1]
        col_permutations = np.argsort(np.random.rand(*arr.shape), axis=0)

        # Use advanced indexing to reorder elements within each column
        # np.arange(N_cols)[None, :] creates a shape (1, N_cols) array for correct indexing
        arr_fully_permuted = arr_permuted_cols[col_permutations, np.arange(N_cols)[None, :]]
        # Create a new DataFrame using the original index and columns
        df_final = pd.DataFrame(
            arr_fully_permuted, 
            index=cdt_sorted.index, 
            columns=cdt_sorted.columns
        )
        permut_cnv_score =CNVScore._cal_cnv_score(df1_bulk, df_final,row_num_np)
        return permut_cnv_score

    @staticmethod
    def permut_scores(pmtimes,cdt_sorted,df1_bulk,row_num_np):
        """Calculate permutation CNV scores"""
        permuted_scores = []
        for i in range(pmtimes):
            permut_cnv_score = CNVScore.get_single_permut_cnv_score(cdt_sorted,df1_bulk,row_num_np)
            permuted_scores.append(permut_cnv_score)
        df_arr = np.array(permuted_scores)
        permuted_scores_df = pd.DataFrame(df_arr.T)
        return permuted_scores_df
