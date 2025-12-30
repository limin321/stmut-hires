import os
import pandas as pd

class BarcodesSortedProcessor:
    """Responsible for processing sorted barcodes"""

    @staticmethod
    def process_barcode_group(sorted_barcodes,parquet_file,ensembl, output_dir, cutoff):
        """Merge spatially adjacent barcodes that collectively meet expression cutoff"""
        merged_barcodes = []
        actual_non_zero = 0

        # Add barcodes one by one and check merged expression
        for barcode in sorted_barcodes:
            merged_barcodes.append(barcode)

            # Read expression matrix for current set of barcodes
            merged_expr = parquet_file.read(columns=merged_barcodes).to_pandas()
            summed_expr = merged_expr.sum(axis=1)  # Combine gene expressions across merged barcodes
            actual_non_zero = int((summed_expr > 0).sum())  # Count genes expressed at least once

            if actual_non_zero >= cutoff:
                break    

        # Warn if below cutoff (only possible when exhausted all barcodes)
        if actual_non_zero < cutoff:
            import warnings
            warnings.warn(
                f"Group {merged_barcodes[0]} reached maximum available barcodes."
                f"(genes: {actual_non_zero}, cutoff: {cutoff})",
                RuntimeWarning
            )
        
        # Save results
        result = pd.concat([ensembl, summed_expr.to_frame(name='merged')], axis=1)
        result.to_csv(
            os.path.join(output_dir,f"{merged_barcodes[0]}.txt"),
            sep="\t", index=False, header=False
        )
        return merged_barcodes, actual_non_zero

