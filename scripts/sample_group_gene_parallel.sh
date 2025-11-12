#! /usr/bin/bash

source activate scvi-env

indir="xx/Visium/stmut_python/BD17_bin8_inputs"
outdir="xx/Visium/stmut_python/two_steps_out"

python ../main.py \
    --exp_h5 ${indir}/filtered_feature_bc_matrix.h5 \
    --cluster_file ${indir}/Graph-Based.csv \
    --spatial_file ${indir}/spatial/tissue_positions.parquet \
    --output_dir ${outdir} \
    --cutoff 3000 \
    --num_processes 3



