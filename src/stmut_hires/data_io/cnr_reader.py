import pandas as pd


class CnrReader:
    """Reads and filters CNR files."""

    def __init__(self, cnr):
        self.cnr = cnr

    def read_cnr(self):
        """Read a CNR file, dropping MT and Y chromosomes."""
        cnr = pd.read_table(self.cnr, sep="\t")
        cnr_name = self.cnr.split('/')[-1]
        cnr = cnr[~cnr['chromosome'].isin(['MT', 'Y'])].copy()
        cnr.loc[cnr['chromosome'] == 'X', 'chromosome'] = 23
        return cnr_name, cnr
