# stmut-hires

## Description
A command-line tool for CNV analysis in high-resolution spatial transcriptomics. It is the upgraded version of [stmut](https://github.com/limin321/stmut) to handle high resolution spatial transcriptomics data (VisiumHD, StereoSeq)

## Getting Started
### Performance Benchmarks
The following benchmarks were recorded using a high-resolution (~680k) barcodes VisiumHD dataset.

|Metric              | Full Pipeline Run (run)                     |
|--------------------|---------------------------------------------|
|Dataset Type        |Visium HD                                    |
|Dataset Size        |682,346 barcodes(Tissue-filtered)            |
|Gene Count          |18,085 Genes                                 |
|Wall Time(Duration) |15h 56m                                      |
|Peak Memory(RAM)    |44.58 GB                                     |
|Validated Hardware  |Intel(R) Xeon(R) Gold 5222 (16 Logical CPUs) |

Memory Scaling: It took ~45GB for 682k spots. This suggests approximately 65–70 MB of RAM per 1,000 spots.  
Recommendation: For a full-resolution 11-million bin dataset, users will need to utilize the binning features (8µm or 16µm) or request high-memory nodes (512GB+) to avoid OOM (Out of Memory) crashes.



### 🖥️ Validated Testing Environment
To ensure reproducibility, performance and stability have been validated on the following enterprise-grade configuration:
- - Hardware Specifications
* **CPU:** Intel(R) Xeon(R) Gold 5222 @ 3.80GHz (Dual-Socket)
* **Architecture:** x86_64 (16 logical CPUs, 8 physical cores)
* **Cache/Tech:** 16.5MB L3 Cache | VT-x Virtualization enabled
* **Memory Management:** NUMA-aware (2 nodes) optimized for parallel processing
- - Operating System & Software
* **OS:** CentOS Linux 7 (Core)
* **Kernel:** x86_64 Baseline
* **Glibc Version:** 2.17 (Manylinux_2_17 compatible)

#### Prerequisites
* Python: 3.9+    
* Key Dependencies: pandas(v2.3.2), cnvkit(v0.9.12), pomegranate(v0.14.8), biopython(v1.86), h5py(v3.15.1)  

### Installation
1.  Clone the repository and create the environment:
```
git clone https://github.com/limin321/stmut-hires.git
cd stmut-hires
conda env create -f environment.yml
conda activate stmut-hires
```
2. Fix Library Link (Required for HPC stability)
```
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
chmod +x $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
```

3. Verify installation:
```
conda activate stmut-hires
stmut-hires --help
```
If successful, you will see the CNV Analysis Pipeline help message.

#### Docker Deployment
We also provided a Dockerfile. You can create your own docker image while inside `stmut-hires` directory.
`docker build -t stmut-hires:v0.1.0 .`
or pull the docker image from Docker Hub [stmut-hires](https://hub.docker.com/r/limin321/stmut-hires) by
`docker pull limin321/stmut-hires:v0.1.0`

Here is an example of running  stmut-hires using docker container:
```
#! /bin/bash
set -e

indir="xx/xxbin8_inputs"
outdir="xx/xx/outs1"
mkdir -p ${outdir}
sudo chmod -R 777 ${outdir}

docker run --rm \
    -v ${indir}:/data/input \
    -v ${outdir}:/data/outs \
    stmut-hires:v0.1.0 run \
    --exp_h5 /data/input/filtered_feature_bc_matrix.h5 \
    --cluster_file /data/input/Graph-Based.csv \
    --spatial_file /data/input/spatial/tissue_positions.parquet \
    --output_dir /data/outs \
    --manul_cutoff 50 \
    --cutoff 3000 \
    --num_processes 10 \
    --annotate_file /data/input/annotate.csv \
    --bulkCNV_file /data/input/bulkCNV.csv \
    --cores 10 \
    --pmtimes 50 \
    --ncluster 6


echo "Finished"
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

annotate.csv -- A two-column tumor/normal annotation csv file. (Case sensitive, 'cluster' NOT 'Cluster' or 'CLUSTER')
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
# Define paths
INDIR="xx/xx/inputs"
OUTDIR="xx/xx/outs"

stmut-hires run \
    --exp_h5 ${INDIR}/filtered_feature_bc_matrix.h5 \
    --cluster_file ${INDIR}/Graph-Based.csv \
    --spatial_file ${INDIR}/spatial/tissue_positions.parquet \
    --output_dir ${OUTDIR} \
    --cutoff 3000 \
    --num_processes 10 \
    --annotate_file ${INDIR}/annotate.csv \
    --bulkCNV_file ${INDIR}/bulkCNV.csv \
    --manul_cutoff 50 \
    --cores 10 \
    --pmtimes 50 \
    --ncluster 6

```

#### Expected output structure
A successful run will include all the following folders. There are too many files in txt, cnr, wtcnr folder, which are not shown here for better demonstration.
```
(base) [stmut_hires_test_proj]$ tree -L 2 outs
outs
├── cdt
│   └── grpWt.cdt
├── cluster_exp
│   ├── Cluster1.parquet
│   ├── Cluster2.parquet
│   ├── Cluster3.parquet
│   ├── Cluster4.parquet
│   ├── Cluster5.parquet
│   ├── Cluster6.parquet
│   ├── Cluster7.parquet
│   ├── Cluster8.parquet
│   ├── Cluster9.parquet
│   └── ensembl.csv
├── cluster_summary
│   ├── Cluster1_barcode_grouping_info.csv
│   ├── Cluster2_barcode_grouping_info.csv
│   ├── Cluster3_barcode_grouping_info.csv
│   ├── Cluster4_barcode_grouping_info.csv
│   ├── Cluster5_barcode_grouping_info.csv
│   ├── Cluster6_barcode_grouping_info.csv
│   ├── Cluster7_barcode_grouping_info.csv
│   ├── Cluster8_barcode_grouping_info.csv
│   └── Cluster9_barcode_grouping_info.csv
├── cnr
├── figures
│   ├── barcodes_counts_histogram.pdf
│   ├── CNVs_OrganizedByGEcluster_UMIcount.pdf
│   ├── CNVs_RankedBySimilarityToDNA_CNVscoreHistogram.pdf
│   ├── CNVs_RankedBySimilarityToDNA_QQplot.pdf
│   └── unrooted_CNVs_clustered_heatmap_class_6clusters.pdf
├── tables
|   ├── caseCNVScore.parquet
|   ├── cluster_barcodes_summary.csv
|   ├── CNVs_OrganizedByGEcluster_UMIcount.cdt
|   ├── CNVs_RankedBySimilarityToDNA.cdt.parquet
|   ├── CNVs_RankedBySimilarityToDNA_CNVscoreHistogram.csv
|   ├── CNVs_RankedbySimilaritytoDNA_Quintiles4Loupe.csv
|   └── permut_CNVscores.parquet
├── txt
└── wtcnr
```

#### 3. Rerun CNV calling step
It is typically difficult to infer copy number alterations on the X-chr with gene expression data. One copy of the X-chr is silenced via X-inactivation. If the inactive copy of X is subjected to a CNA, it would not show up in the gene expression data. If the active copy is gained, better to include X-chr just to rerun the RankedBySimilarity analysis.

Assuming you have run the full pipeline, just need to modify your "bulkCNV.csv" file, and rerun the following script.
**Note**, the outputs inside the `figures` and `tables` directory from full pipeline will overwritten. Please make your own copy if you want to keep the full pipeline outputs.

```
# Define paths
INDIR="xx/xx/inputs"
OUTDIR="xx/xx/outs"

stmut-hires call-cnv \
    --cluster_file ${INDIR}/Graph-Based.csv \
    --output_dir ${OUTDIR} \
    --annotate_file ${INDIR}/annotate.csv \
    --bulkCNV_file ${INDIR}/bulkCNV.csv \
    --pmtimes 50 \
    --ncluster 6
```

#### 4. Clean output
After running the pipeline, you can clean the output by running `stmut-hires clean` subcommand.
The only input is the `--output_dir` you used to run the pipeline.
```
$ stmut-hires clean --help
usage: stmut-hires clean [-h] --output_dir OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Output directory to clean.
```



### Help
```
usage: stmut-hires [-h] {run,call-cnv,clean} ...

CNV Analysis Pipeline.

positional arguments:
  {run,call-cnv,clean}  Available commands
    run                 Run full pipeline.
    call-cnv            Run Step 6: calling CNV..only
    clean               Clean intermediate output folders, keeping only figures/ and tables/.

options:
  -h, --help            show this help message and exit
(stmut-filter) [lchen@localhost stmut-hires]$ stmut-hires run --help
usage: stmut-hires run [-h] --output_dir OUTPUT_DIR --cluster_file CLUSTER_FILE [--dry_run] --exp_h5 EXP_H5 --spatial_file
                       SPATIAL_FILE [--manual_cutoff MANUAL_CUTOFF] [--num_processes NUM_PROCESSES] [--cutoff CUTOFF]
                       [--window WINDOW] [--cores CORES] [--bw_method BW_METHOD] [--annotate_file ANNOTATE_FILE]
                       [--bulkCNV_file BULKCNV_FILE] [--pmtimes PMTIMES] [--ncluster NCLUSTER] [--distance_metric DISTANCE_METRIC]
                       [--linkage_method LINKAGE_METHOD]

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
  --manual_cutoff MANUAL_CUTOFF
                        Cutoff to filter-out barcodes with low gene counts (INT, default: None). You either set this parameter or
                        `--bw_method` to filter-out low-quality barcodes.
  --num_processes NUM_PROCESSES
                        The number of processes used for parallel (default: None)
  --cutoff CUTOFF       The min number of genes in a new spot as grouping cutoff. (default: 1000)
  --window WINDOW       The nearest-neighbor spots for selecting grouping candidates. (default): 100
  --cores CORES         Number of cores to run weighted-median parallelly
  --bw_method BW_METHOD
                        Bandwidth scaling factor for KDE valley detection (e.g., 0.1 to 0.9). For it to work, `--manual_cutoff`
                        needs to be default value None.
  --annotate_file ANNOTATE_FILE
                        Two-column clusters annotated csv file.
  --bulkCNV_file BULKCNV_FILE
                        [Optional] Two-column csv storing bulk-CNV info.
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

### License
This project is licensed under the **MIT License**.  
See the [LICENSE](https://github.com/limin321/stmut-hires/blob/master/LICENSE) file for the full license text.

### Citation
Chen, L., Chang, D., Tandukar, B. et al. STmut: a framework for visualizing somatic alterations in spatial transcriptomics data of cancer. Genome Biol 24, 273 (2023). https://doi.org/10.1186/s13059-023-03121-6


