import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import os
import pyarrow.parquet as pq
from tqdm.auto import tqdm # For automatic selection of notebook/console bar
import time # Import the time module

#import dask.dataframe as dd
#from dask.diagnostics import ProgressBar

def merge_barcodes(spatial_df, expression_file, ensembl_file, output_dir, cutoff=1000, window=100):
    start_time_total = time.time() # Start overall timer
    

    """Process large expression files efficiently with memory optimization """
    cluster = os.path.splitext(os.path.basename(expression_file))[0]
    ensembl = pd.read_csv(ensembl_file) # read ensembl file
    spatial_df = spatial_df.copy() # Load spatial data and build NN index once - to avoid building every grouping
    parquet_file = pq.ParquetFile(expression_file) # only load the meta data without entire data.
    expression_stats, summed_expr = expr_summary(spatial_df,parquet_file)
    end_time_stats_build = time.time() # End timer for stats build

    initial_spatial_df_len = len(spatial_df)
    tqdm.write(f"Starting main barcode processing for {initial_spatial_df_len} initial spots...")
    start_time_main_loop = time.time() # Start timer for main loop
    pbar = tqdm(total=initial_spatial_df_len, desc="Merging Barcodes") # Initialize tqdm progress bar for the main loop

    # Main processing loop
    merged_groups = []
    group_id =1

    # First check: count the total number of non-zero genes in the expression file of the cluster
    total_non_zero_count = (summed_expr != 0).sum() # Non-zero genes counts
    tqdm.write(f"Total number of non-zero gene counts in {cluster}: {total_non_zero_count}")

    # Case1: The total number of non-zero countes genes in one cluster is less than cutoff. Merge all bins into one cluster.
    if total_non_zero_count <= cutoff:
        new_group, spatial_df = cluster_total_genes_less_cutoff(expression_stats,spatial_df,group_id,ensembl,summed_expr,output_dir)
        merged_groups.append(new_group)

    # Case2: The total number of non-zero countes genes in one cluster is less than cutoff. Grouping is required.
    while len(spatial_df) > 0:
        # Rebuild coordinate arrays and KD-tree with current barcodes
        coords = spatial_df[["array_col", "array_row"]].values
        barcodes = spatial_df["barcode"].values
        top_left = spatial_df.iloc[0] # generate a pandas series
        current_barcode = top_left['barcode']

        if current_barcode not in expression_stats.index:
            print(f"Drop current barcode: {current_barcode}")
            spatial_df = spatial_df.drop(top_left.name) # removes the row from spatial_df that has the index label identified by top_left.name. 
            continue

        # Case2-1: process single barcode. The first bin/spot has non-zero gene count > cutoff:
        current_count = expression_stats.loc[current_barcode, 'non_zero_count']
        if current_count >= cutoff:
            process_single_barcode_fast(current_barcode, parquet_file, ensembl, output_dir)
            merged_groups.append({
                'group_id':f"group_{group_id}",
                'main_barcode': current_barcode,
                'total_non_zero':current_count,
                'array_col': top_left['array_col'],
                'array_row': top_left['array_row'],
                'merged_barcodes':current_barcode,
            })

            # remove the grouped bcs
            spatial_df = spatial_df.drop(top_left.name)
            expression_stats = expression_stats.drop(current_barcode, errors='ignore')
            group_id += 1
            pbar.update(1) # Update progress for a single processed barcode
            continue

        # Case2-2: # Need grouping: Find nearest neighbors using pre-build index; 
        current_idx = spatial_df.index.get_loc(top_left.name)
        sorted_barcodes = dynamic_update_window(barcodes, window, cutoff, coords,current_idx,parquet_file)

        # Process barcode group
        merged_barcodes, total_non_zero = process_barcode_group_fast(
            sorted_barcodes,
            expression_stats,
            parquet_file,
            ensembl,
            output_dir,
            cutoff
        )
        if merged_barcodes:
            main_barcode = merged_barcodes[0]
            main_coords = spatial_df.set_index('barcode').loc[main_barcode, ['array_col', 'array_row']]
            
            merged_groups.append({
                'group_id':f"group_{group_id}",
                'main_barcode': main_barcode,
                'total_non_zero': total_non_zero,
                'array_col': main_coords['array_col'],
                'array_row': main_coords['array_row'],
                'merged_barcodes': ','.join(merged_barcodes)
            })

            # Update pbar by the number of barcodes successfully merged into a group
            pbar.update(len(merged_barcodes))

            spatial_df = spatial_df[~spatial_df['barcode'].isin(merged_barcodes)]
            expression_stats = expression_stats.drop(merged_barcodes, errors='ignore')
            group_id += 1
        else:
            spatial_df = spatial_df.drop(top_left.name)
            pbar.update(1) # Update progress for the single dropped item (no group formed)

    pbar.close() # Close the progress bar when the loop completes
    end_time_main_loop = time.time() # End timer for main loop
    print(f"Main barcode processing loop completed in {end_time_main_loop - start_time_main_loop:.2f} seconds.")
   
    # Save grouping information
    pd.DataFrame(merged_groups).to_csv(
        f"{output_dir}/{cluster}_barcode_grouping_info.csv",
        index = False
    )
    end_time_total = time.time() # End overall timer
    print(f"\n--- Script completed in {end_time_total - start_time_total:.2f} seconds total ---")



def expr_summary(spatial_df, parquet_file):
    """ 
    Create expression_stats summary file and summed expression for later grouping.
    
    This function processes the expression data in chunks to be memory efficient.
    It performs two tasks in a single pass:
    1. Builds a summary of non-zero gene counts for each barcode.
    2. Calculates the total expression sum for each gene across the entire dataset.
    
    Args:
        spatial_df (pd.DataFrame): DataFrame with spatial barcode information.
        parquet_file (pyarrow.parquet.ParquetFile): A ParquetFile object pointing 
                                                    to the expression data.

    Returns:
        tuple: A tuple containing:
            - expression_stats (pd.DataFrame): Summary of non-zero counts per barcode.
            - summed_expr (pd.Series): Total expression sum for each gene.
    """
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


def cluster_total_genes_less_cutoff(expression_stats,spatial_df,group_id,ensembl,summed_expr,output_dir):
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
    spatial_df = spatial_df[~spatial_df['barcode'].isin(merged_barcodes)] # this should remove all data making len(spatial_df)=0
    group_id += 1
    # Save results
    result = pd.concat([ensembl, summed_expr.to_frame(name='merged')], axis=1)
    result.to_csv(
        os.path.join(output_dir,f"{merged_barcodes[0]}.txt"),
        sep="\t", index=False, header=False
    )
    return new_group,spatial_df # return this local variable to replace the global variable so it won't execute while loop.


def dynamic_update_window(barcodes, window, cutoff, coords, current_idx, parquet_file):
    """
    Evaluates the sub_total_non_zero_count to dynamically update the window size.
    This guarantees the grouped data meets the cutoff gene count.
    
    The variable `sorted_barcodes` is now initialized before the loop to prevent an UnboundLocalError
    if the while loop condition is not met on the first iteration.
    """
    
    # Initialize nbrs outside the loop
    nbrs = NearestNeighbors(n_neighbors=min(window, len(coords)), algorithm="kd_tree").fit(coords)
    
    # Initial calculation before the loop
    distances, indices = nbrs.kneighbors([coords[current_idx]])
    sorted_indices = indices[0][np.argsort(distances[0])]
    sorted_barcodes = barcodes[sorted_indices]

    # subset expression data and calculate non-zero gene count
    sub_expr = parquet_file.read(columns=sorted_barcodes).to_pandas()
    sub_expr_sum = sub_expr.sum(axis=1)
    sub_total_non_zero_count = (sub_expr_sum != 0).sum()

    # Loop to dynamically update the window and re-evaluate
    while sub_total_non_zero_count < cutoff and window <= len(coords):
        tqdm.write("Dynamically update window size for this cluster.")
        # Dynamic window update logic
        window *= 2  # Double the window size
        tqdm.write(f"New window size is {window}")
        # Re-fit the model with the larger window for the next iteration
        nbrs = NearestNeighbors(n_neighbors=min(window, len(coords)), algorithm="kd_tree").fit(coords)
        
        # Get nearest neighbors with the new window size
        distances, indices = nbrs.kneighbors([coords[current_idx]])
        sorted_indices = indices[0][np.argsort(distances[0])]
        sorted_barcodes = barcodes[sorted_indices]

        # subset expression data and calculate non-zero gene count
        sub_expr = parquet_file.read(columns=sorted_barcodes).to_pandas()
        sub_expr_sum = sub_expr.sum(axis=1)
        sub_total_non_zero_count = (sub_expr_sum != 0).sum()

    return sorted_barcodes



def process_barcode_group_fast(sorted_barcodes, expression_stats, parquet_file, ensembl, output_dir, cutoff):
    """ 
    Merge spatially adjacent barcodes that collectively meet expression cutoff.
    Assumes all sorted_barcodes exist in expression_stats.index
    """

    # Initialize with closest barcode (guaranteed < cutoff from pre-filter)
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

def process_single_barcode_fast(barcode, parquet_file, ensembl, output_dir):
    """Efficiently process a single barcode using columnar reading. """
    # Read just the needed column
    expr = parquet_file.read(columns=[barcode]).to_pandas()
    expr.columns = ['expression']
    result = pd.concat([ensembl, expr], axis=1)
    result.to_csv(
        os.path.join(output_dir, f"{barcode}.txt"),
        sep="\t", index=False, header=False
    )
   

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Merge spatially adjacent barcodes based on gene expression profile.")
    parser.add_argument("--spatial", required=True, help="Spatial coordinate of barcodes, parquet format")
    parser.add_argument("--expression", required=True, help="Spatial gene expression parquet file. No geneID. Columns are barcodes.")
    parser.add_argument("--ensembl", required=True, help="One-column ensembl_ID txt file.")
    parser.add_argument("--output_dir", required=True, help="The output dir to save merged barcodes txt files.")
    parser.add_argument("--cutoff", required=False, type=int, default=1000, help="The min number of genes/Transcripts per merged barcode. (Default: 1000)")
    parser.add_argument("--window", required=False, type=int, default=100, help="The window size to find K nearest neighbor barcodes. (Default: 100)")

    args = parser.parse_args()
    merge_barcodes(
        spatial_df = args.spatial, 
        expression_file = args.expression,
        ensembl_file = args.ensembl, 
        output_dir = args.output_dir, 
        cutoff=args.cutoff, 
        window = args.window
    )

if __name__ == "__main__":
    main()
