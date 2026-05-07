import pandas as pd


class CentromereReader:
    """Reads and processes a centromere BED file."""

    def __init__(self, bed_file):
        self.bed_file = bed_file

    def read_centromere(self):
        centm = pd.read_table(self.bed_file, sep='\t', names=['chromosome', 'start', 'end'])
        centm['chromosome'] = centm['chromosome'].str[3:]
        centm.loc[centm['chromosome'] == 'X', 'chromosome'] = 23
        centm.loc[centm['chromosome'] == 'Y', 'chromosome'] = 24
        centm['chromosome'] = pd.to_numeric(centm['chromosome'])
        centm_sorted = centm.sort_values(by='chromosome', ascending=True)
        centm_sorted['chromosome'] = centm_sorted['chromosome'].astype(str)
        return centm_sorted
