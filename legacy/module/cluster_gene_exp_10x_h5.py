#! /usr/bin/env python

import os
import pandas as pd
import h5py
import numpy as np
import scipy.sparse
from tqdm import tqdm

def cluster_exp_from_10x_h5(h5file, clusterf, output_dir):
    output_dir = os.path.join(output_dir, "cluster_exp")
    os.makedirs(output_dir, exist_ok=True)

    # Load cluster information
    clusterData = pd.read_csv(clusterf)
    clusterData["Graph-based"] = clusterData["Graph-based"].str.replace(" ", "")
    cluster_dict = {
        cl: clusterData.loc[clusterData["Graph-based"] == cl, "Barcode"].tolist()
        for cl in clusterData["Graph-based"].unique()
    }

    # Open HDF5 expression file
    with h5py.File(h5file, 'r') as f:
        matrix = f['matrix']

        # Read matrix structure
        data = matrix['data'][:]
        indices = matrix['indices'][:]
        indptr = matrix['indptr'][:]
        shape = matrix['shape'][:]
        barcodes = matrix['barcodes'][:].astype(str)

        # Try reading gene names
        if 'name' in matrix['features']:
            genes = matrix['features']['id'][:].astype(str)
        else:
            raise ValueError("Can't find gene names in 'features'.")

        # Save gene IDs once
        pd.DataFrame({'gene_id': genes}).to_csv(
            os.path.join(output_dir, "ensembl.csv"),
            index=False
        )

        # Build sparse matrix
        expr = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)

        # Barcode index lookup
        barcode_index_map = {bc: i for i, bc in enumerate(barcodes)}

        # Process each cluster
        for cluster, cluster_barcodes in tqdm(cluster_dict.items(), desc="Processing clusters"):
            valid_indices = [barcode_index_map[bc] for bc in cluster_barcodes if bc in barcode_index_map]
            if not valid_indices:
                print(f"Cluster {cluster}: No valid barcodes.")
                continue

            sub_expr = expr[:, valid_indices].toarray()  # shape: [genes x barcodes]
            sub_barcodes = [barcodes[i] for i in valid_indices]

            df = pd.DataFrame(sub_expr, columns=sub_barcodes)
            # df.insert(0, "gene_id", genes) # Don't save gene_id column to expression file
            #df.to_csv(os.path.join(output_dir, f"{cluster}.csv"), index=False)

            df.to_parquet(
                os.path.join(output_dir, f"{cluster}.parquet"),
                index=False,
                engine='pyarrow',
                compression='snappy'
            )

    print("Finished processing all clusters.")

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate gene expression file for each cluster.")
    parser.add_argument("--exp_h5", required=True, help="Gene expression h5 file. For 10X, filtered_feature_bc_matrix.h5")
    parser.add_argument("--cluster_file", required=True, help="Spatial barcodes Clusters csv file.")
    parser.add_argument("--output_dir", required=True, help="The dir path where the cluster gene expression csv file is saved to.")
    args = parser.parse_args()

    cluster_exp_from_10x_h5(
        h5file=args.exp_h5, 
        clusterf=args.clusterf, 
        output_dir=args.output_dir
        )

if __name__ == "__main__":
    main()

"""
file_dir = "/Volumes/shainlab/Limin/visiumHD_stmut/BD_17/BD_17_bin64outs/outs/binned_outputs/square_008um"
clusterf = os.path.join(file_dir, "Graph-Based.csv")
output_dir="/Volumes/shainlab/Limin/visiumHD_stmut/BD_17/stmut_python/dask_tmp"
h5file = os.path.join(file_dir, "filtered_feature_bc_matrix.h5")

cluster_exp_from_10x_h5(h5file, clusterf, output_dir)

"""
