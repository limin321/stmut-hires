import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import os
#import dask.dataframe as dd
#from dask.diagnostics import ProgressBar

def merge_barcodes(spatial_df, expression_file, ensembl_file, output_dir, cutoff=2000, window=100):
    """Process large expression files in chunks to conserve memory."""
    os.makedirs(output_dir, exist_ok=True)
    ensembl = pd.read_csv(ensembl_file)
    
    # Initialize variables
    merged_groups = []
    group_id = 1
    
    # Process expression data in chunks
    clst=expression_file.strip().split("/")
    cluster = clst[len(clst)-1].strip().split(".")[0]
    reader=pd.read_parquet(expression_file, engine="pyarrow")
    
    # First pass: Build barcode index and expression summary
    barcode_stats = {}

    for barcode in reader.columns[:]:  # Skip gene_id column
        current_count = (reader[barcode] > 0).sum()
        if barcode in barcode_stats:
            barcode_stats[barcode] += current_count
        else:
            barcode_stats[barcode] = current_count
    
    # Convert to DataFrame for easier manipulation
    expression_stats = pd.DataFrame.from_dict(barcode_stats, orient='index', columns=['non_zero_count'])
    
    # Main processing loop
    while len(spatial_df) > 0:
        top_left = spatial_df.sort_values(['array_row', 'array_col']).iloc[0]
        current_barcode = top_left["barcode"]
        
        if current_barcode not in expression_stats.index:
            spatial_df = spatial_df[spatial_df['barcode'] != current_barcode]
            continue
            
        current_count = expression_stats.loc[current_barcode, 'non_zero_count']
        
        if current_count >= cutoff:
            # Process single barcode
            process_single_barcode(current_barcode, expression_file, ensembl, output_dir)
            
            merged_groups.append({
                'group_id': f"group_{group_id}",
                'main_barcode': current_barcode,
                'merged_barcodes': current_barcode,
                'array_col': top_left['array_col'],
                'array_row': top_left['array_row'],
                'total_non_zero': current_count
            })
            
            spatial_df = spatial_df[spatial_df['barcode'] != current_barcode]
            expression_stats = expression_stats.drop(current_barcode)
            group_id += 1
            continue
        
        # Find nearest neighbors
        coords = spatial_df[['array_col', 'array_row']].values
        barcodes = spatial_df['barcode'].values
        n_neighbors = min(window, len(spatial_df))
        nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(coords)
        
        current_idx = np.where(barcodes == current_barcode)[0][0]
        distances, indices = nbrs.kneighbors([coords[current_idx]])
        sorted_indices = indices[0][np.argsort(distances[0])]
        sorted_barcodes = barcodes[sorted_indices]
        
        # Process barcode group
        merged_barcodes, total_non_zero = process_barcode_group(
            sorted_barcodes, 
            expression_stats,
            expression_file,
            ensembl,
            output_dir,
            cutoff
        )
        
        if merged_barcodes:
            main_barcode = merged_barcodes[0]
            main_coords = spatial_df[spatial_df['barcode'] == main_barcode][['array_col', 'array_row']].iloc[0]
            
            merged_groups.append({
                'group_id': f"group_{group_id}",
                'main_barcode': main_barcode,
                'merged_barcodes': ','.join(merged_barcodes),
                'array_col': main_coords['array_col'],
                'array_row': main_coords['array_row'],
                'total_non_zero': total_non_zero
            })
            
            spatial_df = spatial_df[~spatial_df['barcode'].isin(merged_barcodes)]
            expression_stats = expression_stats.drop(merged_barcodes, errors='ignore')
            group_id += 1
    
    # Save grouping information
    pd.DataFrame(merged_groups).to_csv(
        f"{output_dir}/{cluster}_barcode_grouping_info.csv",
        index=False
    )

def process_single_barcode(barcode, expression_file, ensembl, output_dir):
    """Process a single barcode's expression data. No grouping is performed in this case."""
    expr = None
    chunk_parquet = pd.read_parquet(expression_file, engine="pyarrow")
    if barcode in chunk_parquet.columns:
        if expr is None:
            expr = chunk_parquet[[barcode]].copy()
        else:
            expr = pd.concat([expr, chunk_parquet[[barcode]]])
    
    if expr is not None:
        expr.columns = ['expression']
        expr = pd.concat([ensembl, expr], axis=1)
        expr.to_csv(
            os.path.join(output_dir, f"{barcode}.txt"),
            sep='\t', index=False, header=False
        )

def process_barcode_group(barcodes, expression_stats, expression_file, ensembl, output_dir, cutoff):
    """Process a group of barcodes to merge their expression."""
    merged_expr = None
    merged_barcodes = []
    
    for barcode in barcodes:
        if barcode not in expression_stats.index:
            continue
            
        # Load this barcode's expression
        expr = None
        chunk_parquet = pd.read_parquet(expression_file, engine="pyarrow")
        if barcode in chunk_parquet.columns:
            if expr is None:
                expr = chunk_parquet[[barcode]].copy()
            else:
                expr = pd.concat([expr, chunk_parquet[[barcode]]])
        
        if expr is None:
            continue
            
        # Initialize or merge expression
        if merged_expr is None:
            merged_expr = expr.copy()
            merged_expr.columns = ['merged']
        else:
            merged_expr['merged'] = merged_expr['merged'].add(expr.iloc[:, 0], fill_value=0)
        
        merged_barcodes.append(barcode)
        current_non_zero = (merged_expr['merged'] > 0).sum()
        #print(f"The current_non_zero has counts: {current_non_zero}")
        
        if current_non_zero >= cutoff:
            main_barcode = merged_barcodes[0]
            # Add gene_ids back before saving
            result = pd.concat([ensembl, merged_expr], axis=1)
            result.to_csv(
                os.path.join(output_dir, f"{main_barcode}.txt"),
                sep='\t', index=False, header=False
            )
            return merged_barcodes, current_non_zero 
        
    if len(merged_barcodes) > 0 and current_non_zero < cutoff and len(merged_barcodes) == len(barcodes):
        main_barcode = merged_barcodes[0]
        result = pd.concat([ensembl, merged_expr], axis=1)
        result.to_csv(
            os.path.join(output_dir, f"{main_barcode}.txt"),
            sep='\t', index=False, header=False
        )
        return merged_barcodes, current_non_zero
    return [], 0

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
    )

if __name__ == "__main__":
    main()
