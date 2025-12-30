import pandas as pd

class IntegrateAnnotate:
    """Integrate tumor/normal for each barcode"""
    def __init__(self, annotate, merged_barcode_summary):
        self.annotate = annotate
        self.merged_barcode_summary = merged_barcode_summary

    def annotate_barcode(self):
        merged_df = pd.merge(self.merged_barcode_summary, self.annotate, on='cluster', how='left')
        annotated_barcode = merged_df[['barcode', 'cluster', 'annotate','TotalRDs']]
        return annotated_barcode
