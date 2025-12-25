import pandas as pd
class LoadAnnotate:
    """Responsible for cluster annotation load"""
    def load_annotate(annotate_file):
        ann = pd.read_csv(annotate_file)
        ann['cluster'] = ann['cluster'].str.replace(' ', '')
        return ann
