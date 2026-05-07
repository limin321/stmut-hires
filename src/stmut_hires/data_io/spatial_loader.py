import pandas as pd


class SpatialDataLoader:
    """Loads spatial barcode coordinate data from a parquet file."""

    def __init__(self, spatial_file):
        self.spatial_file = spatial_file
        self.coords = None

    def load_spatial_data(self):
        """Load spatial coordinates data."""
        self.coords = pd.read_parquet(self.spatial_file, columns=['barcode', 'array_row', 'array_col'])
        return self.coords.copy()
