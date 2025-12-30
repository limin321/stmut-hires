import pandas as pd 

# 1. Prepare chromosome side-bar
class SortedChrom:
    """Prepare chromosome data for side-bar plot"""
    @staticmethod
    def gene_count_per_chr(cdt_meta):
        """Prepare chromosome data for left-side-bar plot"""
        meta = cdt_meta.copy()
        meta[['chr', 'loc', 'gene']] = meta['NAME'].str.split(':', expand=True)
        meta = meta['chr'].value_counts().reset_index(name = 'count')
        meta['chr'] = pd.to_numeric(meta['chr'], errors='coerce')
        chr_meta = meta.sort_values(by='chr', ascending=True)
        chr_meta = chr_meta[::-1]
        return chr_meta
