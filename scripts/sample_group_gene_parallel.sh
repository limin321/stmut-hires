#! /usr/bin/bash

source activate scvi-env
# python ./scripts/group_barcodes.py \
#     --expression_file ./BD17_bin8_inputs/filtered_feature_bc_matrix.h5 \
#     --cluster_file ./BD17_bin8_inputs/Graph-Based.csv \
#     --spatial_file ./BD17_bin8_inputs/spatial/tissue_positions.parquet \
#     --output_dir /stomics_data/liminData/Visium/stmut_python/BD17_bin8_inputs/


time python ../python/group_barcodes_multiprocess.py \
    --expression_file ./BD17_bin8_inputs/filtered_feature_bc_matrix.h5 \
    --cluster_file ./BD17_bin8_inputs/Graph-Based.csv \
    --spatial_file ./BD17_bin8_inputs/spatial/tissue_positions.parquet \
    --output_dir ./BD17_bin8_inputs/parallel4 \
    --cutoff 3000 \
    --num_processes 4


    
