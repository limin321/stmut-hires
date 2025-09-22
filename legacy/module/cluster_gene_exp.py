#! /usr/bin/env python

import os 
import pandas as pd



def cluster_exp(expression, clusterf, output_dir, chunk_size=1000):
    # Extract ensembl_ID
    first_col_file = f"{output_dir}/ensembl.csv"

    with open(first_col_file, 'w') as f_out:
        f_out.write('gene_id\n')

        for chunk in pd.read_csv(expression, usecols=[0], chunksize=chunk_size):
            chunk.to_csv(f_out, header=False, index=False)
    
    # Create cluster_dict
    output_dir = os.path.join(output_dir,"cluster_exp")
    os.makedirs(output_dir, exist_ok=True)
    clusterData = pd.read_csv(clusterf)
    clusterData["Graph-based"] = clusterData["Graph-based"].str.replace(" ","")

    # Create Cluster-Barcodes dict
    cluster_dict = {}
    for cluster in clusterData["Graph-based"].unique():
        df1 = clusterData[clusterData["Graph-based"]==cluster]
        cluster_dict[cluster] = df1["Barcode"]
 
    # subset gene_exp for each cluster
    with open(expression) as f:
        all_columns = next(f).strip().split(',')

    writers = {}
    for cluster, columns in cluster_dict.items():
        valid_columns = [col for col in columns if col in all_columns]

        output_path = f"{output_dir}/{cluster}.csv"
        writers[cluster] = {
            'columns':valid_columns,
            'file': open(output_path, 'w'),
            'first_chunk':True
        }
        print(f"Cluster {cluster} will save {len(valid_columns)} columns")

    for chunk in pd.read_csv(expression, chunksize=chunk_size, low_memory=False):
        for cluster, writer in writers.items():
            cluster_chunk = chunk[writer['columns']]

            cluster_chunk.to_csv(
                writer['file'],
                mode='a',
                header=writer['first_chunk'],
                index=False
            )
            writer['first_chunk']=False

    for writer in writers.values():
        writer['file'].close()

    print("Finished processing all clusters.")

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Divide the gene expression file to each cluster.")
    parser.add_argument("--expression", required=True, help="The gene expression file in csv format.")
    parser.add_argument("--clusterf", required=True, help="Spatial barcodes Cluster csv file.")
    parser.add_argument("--output_dir", required=True, help="The dir path to save cluster gene expression csv file.")
    parser.add_argument("--chunk_size", required=False, type=int, default=1000, help="Chunck_size of pandas reading a file. (Default: 1000)")

    args = parser.parse_args()

    cluster_exp(
        expression=args.expression, 
        clusterf=args.clusterf, 
        output_dir=args.output_dir, 
        chunk_size=args.chunk_size
        )

if __name__ == "__main__":
    main()


"""
file_dir = "/Volumes/shainlab/Limin/visiumHD_stmut/BD_17/BD_17_bin64outs/outs/binned_outputs/square_008um"
clusterf = os.path.join(file_dir, "Graph-Based.csv")

"""

