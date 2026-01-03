# stmut-hires

## Description
A command-line tool for CNV analysis in high-resolution spatial transcriptomicss. It is the upgraded version of [stmut](https://github.com/limin321/stmut) to handle high resolution spatial transcriptomics data (VisiumHD, StereoSeq)

## Getting Started

⚠️ This project is under active development.  
### Prerequisites
* Python: 3.9+    
* Key Dependencies: pandas(v2.3.2), cnvkit(v0.9.12), pomegranate(v0.14.8), biopython(v1.86), h5py(v3.15.1)  
* Tested OS: "CentOS Linux 7 (Core)"   

### Performance Benchmarks
-- add table - software, hardware

### Installation
```
# 1. First create a conda environment
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda config --set channel_priority strict
conda create -n stmut-hires python=3.10 cnvkit pomegranate biopython h5py -c conda-forge -c bioconda

# 2. Activate the new env:
conda activate stmut-hires

# 3.  Clone the stmut-hires to your hpc
git clone https://github.com/limin321/stmut-hires.git
cd stmut-hires
pip install .

# 4. Test install successfully ...
stmut-hires --help
usage: stmut-hires [-h] {run,call-cnv} ...

CNV Analysis Pipeline.

positional arguments:
  {run,call-cnv}  Available commands
    run           Run full pipeline.
    call-cnv      Run Step 6: calling CNV..only

options:
  -h, --help      show this help message and exit
```

### Usage example:
#### 1. Prepare inputs files

|   inputs                      | Note                            |
|-------------------------------|---------------------------------|
| filtered_feature_bc_matrix.h5 | SpaceRanger output              |
| tissue_positions.parquet      | SpaceRanger output              |
| Graph-Based.csv               | Exported from LoupeBrowser      |
| outdir                        | output_dir path                 |
| annotate.csv                  | cluster tumor/normal annotation |
| bulkCNV.csv  (Optional)       | bulk CNV info                   |

annotate.csv -- A two-column tumor/normal annotation csv file.
```
cluster,annotate
cluster 1,tumor
cluster 2,normal
cluster 3,tumor
cluster 4,normal
cluster 5,tumor
```


bulkCNV.csv - A two-column csv file. '1'-gain, '-1'-loss
```
arms,gainloss
1p,1
3p,1
3q,1
5p,-1
```

For other inputs, please refer to the `--help` message

#### 2. Running full pipeline
```
indir="xx/xx/inputs"
outdir="xx/xx/outs"

stmut-hires run \
    --exp_h5 ${indir}/filtered_feature_bc_matrix.h5 \
    --cluster_file ${indir}/Graph-Based.csv \
    --spatial_file ${indir}/spatial/tissue_positions.parquet \
    --output_dir ${outdir} \
    --cutoff 3000 \
    --num_processes 10 \
    --annotate_file ${indir}/annotate.csv \
    --bulkCNV_file ${indir}/bulkCNV.csv \
    --cores 10 \
    --pmtimes 50 \
    --ncluster 6

```

#### 3. Rerun CNV calling step
It is typically difficult to infer copy number alterations on the X-chr with gene expression data. One copy of the X-chr is silenced via X-inactivation. If the inactive copy of X is subjected to a CNA, it would not show up in the gene expression data. If the active copy is gained, better to include X-chr just to rerun the RankedBySimilarity analysis.

Assue you have run the full pipeline, just need to modify your "bulkCNV.csv" file, and rerun the following script.
Note, the outputs inside the `figures` and `tables` directory from full pipeline will overwritten. Please make your own copy if you want to keep the full pipeline outputs.

```
indir="xx/xx/inputs"
outdir="xx/xx/outs"

stmut-hires call-cnv \
    --cluster_file ${indir}/Graph-Based.csv \
    --output_dir ${outdir} \
    --annotate_file ${indir}/annotate.csv \
    --bulkCNV_file ${indir}/bulkCNV.csv \
    --pmtimes 50 \
    --ncluster 6
```

### Help
```
stmut-hires run --help
usage: stmut-hires run [-h] --output_dir OUTPUT_DIR --cluster_file CLUSTER_FILE [--dry_run] --exp_h5 EXP_H5
                       --spatial_file SPATIAL_FILE [--num_processes NUM_PROCESSES] [--cutoff CUTOFF]
                       [--window WINDOW] [--cores CORES] [--annotate_file ANNOTATE_FILE]
                       [--bulkCNV_file BULKCNV_FILE] [--pmtimes PMTIMES] [--ncluster NCLUSTER]
                       [--distance_metric DISTANCE_METRIC] [--linkage_method LINKAGE_METHOD]

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Output directory.
  --cluster_file CLUSTER_FILE
                        Barcodes cluster info.(csv);Graph-based.csv downloaded from Loupe Browser.
  --dry_run             Dry run mode
  --exp_h5 EXP_H5       Gene expression matrix (h5)
  --spatial_file SPATIAL_FILE
                        Barcodes spatial coordinates.(./spatial/tissue_positions.parquet)
  --num_processes NUM_PROCESSES
                        The number of processes used for parallel (default: None)
  --cutoff CUTOFF       The min number of genes in a new spot as grouping cutoff. (default: 1000)
  --window WINDOW       The nearest-neighbor spots for selecting grouping candidates. (default): 100
  --cores CORES         Number of cores to run weighted-median parallelly
  --annotate_file ANNOTATE_FILE
                        Two-column clusters annotated csv file.
  --bulkCNV_file BULKCNV_FILE
                        Two-column csv storing bulk-CNV info.
  --pmtimes PMTIMES     The number of permutation times.
  --ncluster NCLUSTER   Number of clusters for CNV plot.
  --distance_metric DISTANCE_METRIC
                        Distance metric for clustering
  --linkage_method LINKAGE_METHOD
                        Linkage method for hierarchical clustering
```


### Authors
Limin Chen
lynnchen31@gmail.com





