import os
import pandas as pd

class SingleBarcodeProcessor:
    """Responsible for processing single barcodes"""
    @staticmethod
    def process_single_barcode(barcode, parquet_file, ensembl, output_dir):

        """Efficiently process a single barcode using columnar reading"""
        expr = parquet_file.read(columns=[barcode]).to_pandas()
        expr.columns = ['expression']
        result = pd.concat([ensembl, expr], axis=1)
        result.to_csv(
            os.path.join(output_dir, f"{barcode}.txt"),
            sep="\t", index=False, header=False
        )

