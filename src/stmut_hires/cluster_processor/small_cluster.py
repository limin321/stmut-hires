import os
import pandas as pd


class SmallClusterProcessor:
    """Responsible for processing clusters with total genes less than cutoff"""
    @staticmethod
    def process_small_cluster(expression_stats,spatial_df,group_id,ensembl,summed_expr,output_dir):
        """Process cluster where non-zero genes is less than cutoff"""
        merged_barcodes = expression_stats.index.values.tolist()
        actual_non_zero = expression_stats['non_zero_count'].sum()
        top_left = spatial_df.iloc[0]

        new_group={
            'group_id':f"group_{group_id}",
            'main_barcode': merged_barcodes[0],
            'total_non_zero':actual_non_zero,
            'array_col': top_left['array_col'],
            'array_row': top_left['array_row'],
            'merged_barcodes': ','.join(merged_barcodes)
        }

        # Save results
        result = pd.concat([ensembl, summed_expr.to_frame(name='merged')], axis=1)
        result.to_csv(
            os.path.join(output_dir,f"{merged_barcodes[0]}.txt"),
            sep="\t", index=False, header=False
        )
        return new_group

