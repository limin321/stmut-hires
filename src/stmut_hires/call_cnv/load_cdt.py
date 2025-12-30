import pandas as pd

class Loadcdt:
    """Read cdt file"""
    @staticmethod
    def cdt_loader(cdt_file):
        cdt = pd.read_csv(cdt_file, sep="\t", low_memory=False)
        cdt_meta = cdt.iloc[2:,1:3]
        return cdt, cdt_meta
