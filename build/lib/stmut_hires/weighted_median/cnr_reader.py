import pandas as pd

class CnrReader:
    """Responsible for handling cnr files"""
    def __init__(self, cnr):
        self.cnr = cnr

    def read_cnr(self):
        """Read-in cnr file for each spot"""
        cnr = pd.read_table(self.cnr, sep="\t")
        cnr_name = self.cnr.split('/')[-1]
        cnr = cnr[~cnr['chromosome'].isin(['MT', 'Y'])].copy()
        cnr.loc[cnr['chromosome']=='X', 'chromosome'] = 23
        return cnr_name, cnr
