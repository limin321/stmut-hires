import os
import glob
import pandas as pd
import pyarrow.parquet as pq
import multiprocessing

from module.merge_gene_parquet_optimized import merge_barcodes # finish in 11hours
from module.cluster_gene_exp_10x_h5 import cluster_exp_from_10x_h5




def process_cluster(cluster_path, coords, output_dir_txt, ensembl_path, cutoff=1000, window=100):
    print(f"Working on {cluster_path}")
    # Extract all barcodes of the cluster 
    table = pq.read_table(cluster_path)  # loads no rows
    barcodes = table.schema.names

    # Subset the spatial barcodes of the cluster
    cluster_coor = coords[coords['barcode'].isin(barcodes)]

    merge_barcodes(
        spatial_df=cluster_coor, 
        expression_file=cluster_path, 
        ensembl_file=ensembl_path, 
        output_dir=output_dir_txt, 
        cutoff=cutoff, 
        window=window)

def process_cluster_parallel(clusterf, expf, spatialf, output_dir, cutoff=1000, window=100, num_processes=None):
    """
    Proceess clusters in parallel with safe output handling.
    """
    if num_processes is None:
        num_processes = multiprocessing.cpu_count() - 1

    # Step1: genereate expression file for each cluster
    # Generate each cluster_gene_exp file
    cluster_exp_from_10x_h5(expf, clusterf, output_dir)
    print("Each cluster_gene_exp file are generated.")

    # Step2: Group each cluster -- There is a bug, when the input cluster_exp not found, no Error.
    clusters = glob.glob(f"{output_dir}/cluster_exp/Cluster*.parquet")
    coords = pd.read_parquet(spatialf, columns=['barcode', 'array_row', 'array_col'])
    ensembl_path = f"{output_dir}/cluster_exp/ensembl.csv"

    # Create txt folder to save merged barcodes:
    output_dir_txt = f"{output_dir}/txt"
    os.makedirs(output_dir_txt, exist_ok=True)

    # Make sure coords is in the right format (pre-process once)
    coords = coords.copy()

    with multiprocessing.Pool(processes=num_processes) as pool:
        args = [(cluster, coords, output_dir_txt, ensembl_path, cutoff, window) for cluster in clusters]
        pool.starmap(process_cluster, args)

        
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Grouping barcodes within each cluster.")
    parser.add_argument("--expression_file", required=True, help="Gene expression matrix")
    parser.add_argument("--cluster_file", required=True, help="Barcodes cluster info.")
    parser.add_argument("--spatial_file", required=True, help="Barcodes spatial coordinates.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument("--num_processes", required=False, default=None,type=int, help="The number of processes used for parallel")
    parser.add_argument("--cutoff", required=False, default=1000,type=int, help="The number of genes in the new spot used for grouping.")
    parser.add_argument("--window", required=False, default=100,type=int, help="The nearest-neighbor spots for selecting grouping candidates.")

    args = parser.parse_args()


    process_cluster_parallel(
        clusterf=args.cluster_file, 
        expf=args.expression_file,
        spatialf=args.spatial_file, 
        output_dir=args.output_dir, 
        num_processes=args.num_processes,
        cutoff=args.cutoff, 
        window=args.window
    )
    
if __name__ == "__main__":
    main()





"""
file_dir = "/Users/limin/limin_practice/Apps/DNAnexus/stmutCNVtest/scripts/rep2"
clusterf = os.path.join(file_dir, "Graph-Based.csv")
#expf = os.path.join(file_dir,"filtered_feature_bc.csv")
h5file = os.path.join(file_dir, "filtered_feature_bc_matrix.h5")
spatialf = f"{file_dir}/spatial/tissue_positions.csv"
output_dir = "/Users/limin/limin_practice/Apps/DNAnexus/stmutCNVtest/scripts/rep2/stmut_outs/"

dir1="/Volumes/shainlab/Limin/visiumHD_stmut/BD_17/BD_17_bin64outs/outs/binned_outputs/square_008um"
python ../scripts/group_barcodes.py --expression_file ${dir1}/filtered_feature_bc_matrix.h5 --cluster_file ${dir1}/Graph-Based.csv --spatial_file ${dir1}/spatial/tissue_positions.parquet --output_dir /Volumes/shainlab/Limin/visiumHD_stmut/BD_17/stmut_python/dask_tmp

"""