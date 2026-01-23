import pandas as pd
import os
import pyarrow.parquet as pq
from tqdm.auto import tqdm # For automatic selection of notebook/console bar
import time
import logging

# import from current package(processor)
from stmut_hires.cluster_processor.expression_statistic import ExpressionStatisticsCalculator
from stmut_hires.cluster_processor.track_process import ProcessTracker
from stmut_hires.cluster_processor.small_cluster import SmallClusterProcessor
from stmut_hires.cluster_processor.single_barcode import SingleBarcodeProcessor
from stmut_hires.cluster_processor.sort_barcode import BarcodeSorter
from stmut_hires.cluster_processor.merge_sorted_barcode import BarcodesSortedProcessor

class BarcodeMerger:
    """Main class responsible for coordinating the barcode merging process"""
    
    def __init__(self, cluster_spatial_coords, expression_file, ensembl_file, output_dir, cutoff=1000, window=100):
        self.cluster_spatial_coords = cluster_spatial_coords
        self.expression_file = expression_file
        self.ensembl_file = ensembl_file
        self.output_dir = output_dir
        self.cutoff = cutoff
        self.window = window
        self.logger = logging.getLogger(__name__)
        
    def merge_barcodes(self):
        """Main method to merge barcodes"""
        start_time_total = time.time()
        
        cluster = os.path.splitext(os.path.basename(self.expression_file))[0]
        ensembl = pd.read_csv(self.ensembl_file)
        cluster_spatial_coords = self.cluster_spatial_coords.copy()
        parquet_file = pq.ParquetFile(self.expression_file)
        
        # Calculate expression statistics
        stats_calculator = ExpressionStatisticsCalculator()
        expression_stats, summed_expr = stats_calculator.calculate_expression_stats(cluster_spatial_coords, parquet_file)
        
        initial_spatial_df_len = len(cluster_spatial_coords)
        
        tqdm.write(f"Starting main barcode processing for {initial_spatial_df_len} initial spots...")
        start_time_main_loop = time.time()
        
        # Initialize progress tracker
        progress_tracker = ProcessTracker(initial_spatial_df_len, "Merging Barcodes")
        
        # Main processing
        merged_groups = []
        group_id = 1
        
        total_non_zero_count = (summed_expr != 0).sum()
        tqdm.write(f"Total number of non-zero gene counts in {cluster}: {total_non_zero_count}")

        # Case 1: Total non-zero genes <= cutoff
        if total_non_zero_count <= self.cutoff:
            small_processor = SmallClusterProcessor()
            new_group = small_processor.process_small_cluster(
                expression_stats, cluster_spatial_coords, group_id, ensembl, summed_expr, self.output_dir
            )
            merged_groups.append(new_group)
            cluster_spatial_coords = cluster_spatial_coords[~cluster_spatial_coords['barcode'].isin(expression_stats.index.values.tolist())]
            group_id += 1

        # Case 2: Process remaining barcodes
        while len(cluster_spatial_coords) > 0:
            result = self._process_next_barcode(
                cluster_spatial_coords, expression_stats, parquet_file, ensembl, 
                group_id, progress_tracker
            )
            
            if result:
                merged_group, updated_spatial_df, updated_expression_stats, updated_group_id = result
                merged_groups.append(merged_group)
                cluster_spatial_coords = updated_spatial_df
                expression_stats = updated_expression_stats
                group_id = updated_group_id

        progress_tracker.close()
        end_time_main_loop = time.time()
        self.logger.info(f"Main barcode processing loop completed in {end_time_main_loop - start_time_main_loop:.2f} seconds.")
        
        # Save grouping information
        self._save_grouping_info(merged_groups, cluster)
        
        end_time_total = time.time()
        self.logger.info(f"\n--- Script completed in {end_time_total - start_time_total:.2f} seconds total ---")
    
    def _process_next_barcode(self, cluster_spatial_coords, expression_stats, parquet_file, ensembl, group_id, progress_tracker):
        """Process the next barcode in the spatial dataframe"""
        top_left = cluster_spatial_coords.iloc[0]
        current_barcode = top_left['barcode']

        if current_barcode not in expression_stats.index:
            self.logger.info(f"Drop current barcode: {current_barcode}")
            cluster_spatial_coords = cluster_spatial_coords.drop(top_left.name)
            progress_tracker.update(1)
            return None, cluster_spatial_coords, expression_stats, group_id

        current_count = expression_stats.loc[current_barcode, 'non_zero_count']
        # Case 2-1: Single barcode meets cutoff
        if current_count >= self.cutoff:
            single_processor = SingleBarcodeProcessor()
            single_processor.process_single_barcode(current_barcode, parquet_file, ensembl, self.output_dir)
            
            merged_group = ({
                'group_id': f"group_{group_id}",
                'main_barcode': current_barcode,
                'total_non_zero': current_count,
                'array_col': top_left['array_col'],
                'array_row': top_left['array_row'],
                'merged_barcodes': current_barcode,
            })

            cluster_spatial_coords = cluster_spatial_coords.drop(top_left.name)
            expression_stats = expression_stats.drop(current_barcode, errors='ignore')
            progress_tracker.update(1)
            group_id += 1
            return merged_group, cluster_spatial_coords, expression_stats, group_id

        # Case 2-2: Need grouping
        coords = cluster_spatial_coords[["array_col", "array_row"]].values
        barcodes = cluster_spatial_coords["barcode"].values
        current_idx = cluster_spatial_coords.index.get_loc(top_left.name)
        
        grouper = BarcodeSorter()
        sorted_barcodes = grouper.dynamic_update_window(
            barcodes, self.window, self.cutoff, coords, current_idx, parquet_file
        )

        group_processor = BarcodesSortedProcessor()
        merged_barcodes, total_non_zero = group_processor.process_barcode_group(
            sorted_barcodes, parquet_file, ensembl, self.output_dir, self.cutoff
        )
        
        if merged_barcodes:
            main_barcode = merged_barcodes[0]
            main_coords = cluster_spatial_coords.set_index('barcode').loc[main_barcode, ['array_col', 'array_row']]
            
            merged_group = {
                'group_id': f"group_{group_id}",
                'main_barcode': main_barcode,
                'total_non_zero': total_non_zero,
                'array_col': main_coords['array_col'],
                'array_row': main_coords['array_row'],
                'merged_barcodes': ','.join(merged_barcodes)
            }

            progress_tracker.update(len(merged_barcodes))
            cluster_spatial_coords = cluster_spatial_coords[~cluster_spatial_coords['barcode'].isin(merged_barcodes)]
            expression_stats = expression_stats.drop(merged_barcodes, errors='ignore')
            group_id += 1
            return merged_group, cluster_spatial_coords, expression_stats, group_id

    
    def _save_grouping_info(self, merged_groups, cluster):
        """Save grouping information to CSV"""
        pd.DataFrame(merged_groups).to_csv(
            f"{self.output_dir}/{cluster}_barcode_grouping_info.csv",
            index=False
        )

