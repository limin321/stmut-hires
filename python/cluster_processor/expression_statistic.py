import pandas as pd

class ExpressionStatisticsCalculator:
    """Responsible for calculating expression statistics """
    @staticmethod
    def calculate_expression_stats(spatial_df,parquet_file):
        """Create expression_stats summary file and summed_expression for later grouping """
        # Building expression statistics...
        expression_stats = pd.DataFrame(index=spatial_df['barcode'], columns=['non_zero_count'])
        expression_stats['non_zero_count'] = 0

        summed_expr = None

        # Process parquet file in chunks
        for i, batch in enumerate(parquet_file.iter_batches(batch_size=50000)):
            batch_df = batch.to_pandas()

            # Initialize summed_expr with the correct gene index from the first batch
            if i == 0:
                summed_expr = pd.Series(0, index=batch_df.index)

            # Accumulate sums for each gene
            summed_expr += batch_df.sum(axis=1)

            # Accumulate non-zero counts for barcodes that exist in both
            batch_non_zero = (batch_df > 0).sum()
            common_barcodes = batch_non_zero.index.intersection(expression_stats.index)
            expression_stats.loc[common_barcodes, 'non_zero_count'] += batch_non_zero[common_barcodes]

        # Clean up missing barcodes
        expression_stats = expression_stats.dropna()
        return expression_stats, summed_expr

