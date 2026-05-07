import pandas as pd


class AnnotateReader:
    """Reads the cluster annotation CSV file."""

    @staticmethod
    def load_annotate(annotate_file):
        ann = pd.read_csv(annotate_file)
        ann['cluster'] = ann['cluster'].str.replace(' ', '')
        return ann
