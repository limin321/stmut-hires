# Architecture

## Overview

`stmut-hires` is a command-line tool for Copy Number Variation (CNV) analysis on high-resolution spatial transcriptomics data (Visium HD, StereoSeq). It ingests a Cell Ranger gene-expression matrix together with Loupe-Browser cluster assignments and Visium spatial coordinates, then runs a six-stage pipeline that splits expression by cluster, spatially merges neighboring barcodes to meet a per-spot gene-count cutoff, calls bin-level copy-number ratios with CNVkit, smooths them with a weighted median per chromosome arm, and finally classifies barcodes as CNV-positive or negative through permutation against an optional bulk CNV reference.

The package layout (`src/stmut_hires/`) is organized so that each subpackage owns one well-defined pipeline stage or one well-defined responsibility (I/O, plotting, configuration). All subpackage `__init__.py` files are empty; every importer uses the full module path. The dependency graph is strictly tiered — `main.py` orchestrates the workflow modules, the workflow modules call leaf processors, and no two modules import from each other.

## Source layout

```
src/stmut_hires/
├── main.py                    # Pipeline orchestrator and console-script entry point
├── cli_parser.py              # argparse definitions for the run / call-cnv / clean subcommands
├── config.py                  # BaseConfig → InitialStepConfig → MergerConfig → CNVAnalysisConfig
├── data_io/                   # File readers, input validator, output cleaner
├── split_cluster_expr/        # Step 1 — KDE-filter and split expression by cluster
├── cluster_processor/         # Step 2 worker — spatial KD-tree merging per cluster
├── parallel_cluster_runner/   # Step 2 driver — parallelizes cluster_processor across clusters
├── weighted_median/           # Step 4 — per-arm weighted median over cnr files
├── call_cnv/                  # Step 6 — CNV calling, permutation, FDR, quintile assignment
├── visualization/             # All matplotlib/seaborn plots + chromosome-meta helpers
└── data/                      # Bundled hg38 reference files (TSV/BED)
```

## Pipeline

The `run` subcommand executes six steps end-to-end. The `call-cnv` subcommand executes only Step 6 (used to re-run CNV calling after editing `bulkCNV.csv`).

```
Inputs (.h5, Graph-Based.csv, tissue_positions.parquet, annotate.csv, bulkCNV.csv)
    │
    ▼
Step 1  split_cluster_expr.ClusterExpressionProcessor
        Reads the 10x h5 matrix, applies a KDE valley cutoff to drop
        low-gene-count barcodes, splits the surviving expression matrix
        by Loupe cluster.
        → cluster_exp/Cluster*.parquet, cluster_exp/ensembl.csv,
          figures/gene_counts_before_after_filtering.png
    │
    ▼
Step 2  parallel_cluster_runner.ParallelClusterProcessor
        For each cluster parquet, KD-tree-sorts barcodes by spatial
        proximity and progressively merges neighbors until each
        merged group exceeds the gene-count cutoff. One worker per
        cluster (multiprocessing.Pool).
        → txt/<barcode>.txt (one per merged group),
          cluster_summary/<cluster>_barcode_grouping_info.csv
    │
    ▼
Step 3  cnvkit import-rna  (external subprocess)
        Run from txt/ with ulimit -s ulimited; uses bundled
        ensembl-gene-info.hg38.tsv (-g) and
        tcga-skcm.cnv-expr-corr.tsv (-c). main.py invokes cnvlib
        through an inline python -c that monkey-patches
        DataFrame.iteritems → DataFrame.items for pandas-2 support.
        → cnr/*.cnr
    │
    ▼
Step 4  weighted_median.WorkflowOrchestrator
        For each cnr file, joins to the hg38 centromere BED and
        computes a weight-aware median (wquantiles) per chromosome
        arm. Runs in a multiprocessing.Pool sized by --cores.
        → wtcnr/<name>.cnr, wtcnr/summary.csv
    │
    ▼
Step 5  cnvkit.py export cdt  (external subprocess)
        Concatenates wtcnr/*.cnr into a single CDT matrix.
        → cdt/grpWt.cdt
    │
    ▼
Step 6  call_cnv.CNVCallPlotWorkflow
        Loads the CDT, caps log2 ratios to [-1,1], drops rows with
        >5% zeros, subtracts non-tumor cluster medians, merges in
        annotate.csv, sorts by cluster + total reads, draws the
        cluster-sorted heatmap and the unrooted clustered heatmap.
        If bulkCNV.csv is supplied, runs a double-permutation CNV
        score, computes FDR, assigns Q1–Q5 quintiles per barcode.
        → tables/*.csv, tables/*.parquet, figures/*.pdf
```

## Subpackages

### `data_io/`

Pure I/O — no analysis. Every reader is a small class with a single `read_*` or `load_*` method.

| Module | Class | Reads |
|--------|-------|-------|
| `spatial_loader.py` | `SpatialDataLoader` | `tissue_positions.parquet` → `(barcode, array_row, array_col)` |
| `h5_reader.py` | `H5MatrixReader` | `filtered_feature_bc_matrix.h5` → CSC components + barcodes + gene IDs |
| `cluster_reader.py` | `ClusterReader` | `Graph-Based.csv` → `{cluster_name: [barcode, …]}` |
| `cnr_reader.py` | `CnrReader` | one `.cnr` TSV (drops MT/Y, renames X → 23) |
| `centromere_reader.py` | `CentromereReader` | bundled `hg38_centromereSimple.bed` |
| `cdt_loader.py` | `CdtLoader` | `cdt/grpWt.cdt` → `(cdt, cdt_meta)` |
| `annotate_reader.py` | `AnnotateReader` | `annotate.csv` (cluster → tumor/normal) |
| `validators.py` | `InputValidator` | header + existence checks for `run` and `call-cnv` |
| `output_cleaner.py` | `OutputCleaner` | deletes everything under `output_dir` except `figures/` and `tables/` |

### `split_cluster_expr/` (Step 1)

`ClusterExpressionProcessor.process()` is the step-1 driver. It calls `H5MatrixReader` and `ClusterReader`, runs `FilterOutCrappyData.get_filter_cutoff` (a `scipy.stats.gaussian_kde` + `scipy.signal.find_peaks` valley-finder on the per-barcode gene-count distribution), masks the CSC matrix, and hands off to `ExpressionMatrixBuilder` (sparse-matrix helpers) and `GeneDataWriter` (parquet output, snappy-compressed via pyarrow). `QCplot` from `visualization/` produces the before/after gene-count histogram.

### `cluster_processor/` (Step 2 worker)

Per-cluster merging logic. `BarcodeMerger` is the entry point used by the parallel runner; it composes:

- `ExpressionStatisticsCalculator` — streams the parquet in 50 k-row batches to compute per-barcode non-zero counts.
- `BarcodeSorter` — `sklearn.neighbors.KDTree` nearest-neighbor lookup with adaptive window-doubling.
- `BarcodesSortedProcessor` — incrementally merges sorted neighbors until the gene-count cutoff is satisfied; writes one `.txt` per merged group.
- `SingleBarcodeProcessor` — fast path for barcodes that already exceed the cutoff alone.
- `SmallClusterProcessor` — fallback that collapses an entire small cluster into one group.
- `ProcessTracker` — thin tqdm wrapper.

Outputs are written into `txt/`; the grouping-info CSV is moved to `cluster_summary/` after the step completes.

### `parallel_cluster_runner/` (Step 2 driver)

- `OutputDirManager.create_output_dir()` creates `txt/`, `cnr/`, `wtcnr/`, `cdt/`, `figures/`, `cluster_summary/`. `move_csvs_to_summary()` relocates `txt/*.csv` once Step 2 finishes.
- `ParallelClusterProcessor.parallel_cluster_merger()` glob-finds `cluster_exp/Cluster*.parquet` and dispatches `_process_cluster_wrapper` over a `multiprocessing.Pool` (size = `--num_processes` or `cpu_count() - 1`).
- `ClusterMerger.process_cluster(cluster_path)` peeks the parquet schema to recover the cluster's barcodes, subsets the spatial-coords DataFrame, and constructs the `BarcodeMerger` for that cluster.

### `weighted_median/` (Step 4)

- `ArmWeightedMedian.chr_arm_weighted_median()` walks p/q arms using the centromere BED and `wquantiles.median` to compute a weight-aware (log2-ratio × probes) median per arm.
- `WtMedianWriter` writes the per-cnr smoothed file plus the cross-sample `wtcnr/summary.csv`.
- `WorkflowOrchestrator.run_parallel_workflow()` dispatches one worker per cnr file over a `multiprocessing.Pool(config.cores)`.

### `call_cnv/` (Step 6)

`CNVCallPlotWorkflow` is the orchestrator. Its `cnv_workflow(ncluster, bulk_csv, pmtimes, distance_metric, linkage_method)` method drives the full step:

- `ClusterSummary.merge_cluster_summary` concatenates `cluster_summary/Cluster*_barcode_grouping_info.csv` and joins to `Graph-Based.csv`.
- `IntegrateAnnotate.annotate_barcode` injects tumor/normal labels from `annotate.csv`.
- `InferCNV.cdt_processor` caps log2 ratios at ±1, drops rows whose non-tumor zero-fraction exceeds 5 %, and subtracts the non-tumor median gene-wise.
- `InferCNV.sort_cluster_totalreads` orders barcodes by cluster then total reads for the cluster-sorted heatmap.
- `CNVSortedByBulkCNV.run_pipeline` (only if `--bulkCNV_file` is supplied) runs a vectorized double permutation in `CNVScore.permut_scores`, computes empirical FDR in `CNVScoreFDR`, and assigns Q1–Q5 buckets in `AssignBarcodeQuintile`.
- `CNVWriter` (a static-method class) emits every CSV/parquet under `tables/`.

The corresponding plots are produced by `visualization.CNVSortedbyClusters`, `visualization.UnrootedCNVHeatmap`, and `visualization.AnalysisVisualizer`.

### `visualization/`

Every chart the pipeline produces. Each class takes its data in via the constructor and writes a PDF/PNG into `<output_dir>/figures/`:

- `QCplot` — gene-count histogram before/after KDE filtering (Step 1).
- `CNVSortedbyClusters` — `figures/CNVs_OrganizedByGEcluster_UMIcount.pdf`.
- `UnrootedCNVHeatmap` — `seaborn.clustermap` with a chromosome side-bar; supports user-supplied `distance_metric` and `linkage_method`.
- `AnalysisVisualizer` — QQ plot, CNV-score histogram, merged-barcode count histogram (only when bulk CNV is supplied).
- `ClusterTotalReads`, `SortedChrom` — small helpers shared across plotting calls.

### `data/`

Bundled reference assets resolved at runtime via `importlib.resources.files("stmut_hires.data")` in `BaseConfig.__init__`:

- `ensembl-gene-info.hg38.tsv` — gene model passed to `cnvkit import-rna -g` (Step 3).
- `tcga-skcm.cnv-expr-corr.tsv` — TCGA-SKCM CNV-expression correlation table passed to `cnvkit import-rna -c` (Step 3).
- `reference/hg38_centromereSimple.bed` — hg38 centromere coordinates read by `CentromereReader` for Step 4.

## Configuration

A single configuration object flows through the entire pipeline. The hierarchy in `config.py`:

```
BaseConfig
  ├── output_dir, dry_run, cutoff (1000), window (100),
  │   num_processes, cores (4)
  └── gene_info, corr_file, bed_file        ← bundled data paths

InitialStepConfig(BaseConfig)
  ├── clusterf, exp_h5
  ├── bw_method (0.1), filter_cutoff
  └── expression_file, ensembl_file          ← set by set_step1_outputs()

MergerConfig(InitialStepConfig)
  └── spatial_file

CNVAnalysisConfig(MergerConfig)
  ├── annotate_csv, bulk_csv
  ├── pmtimes (5), ncluster (6)
  └── distance_metric ('euclidean'), linkage_method ('ward')
```

`main.py` always constructs a `CNVAnalysisConfig`; CLI arguments not supplied by the active subcommand land as `None` and are simply unused by the modules that don't need them.

## CLI

`cli_parser.CommandLineParser.parse_args()` defines three subcommands sharing argparse parent parsers (`global`, `step1-5`, `step6`):

- `run` — full pipeline. Requires `--output_dir`, `--cluster_file`, `--exp_h5`, `--spatial_file`. Optional flags cover all step-1-5 and step-6 parameters listed in the README's help section.
- `call-cnv` — Step 6 only. Requires `--output_dir`, `--cluster_file`; `--annotate_file` is effectively required at runtime by `InputValidator.validate_call_cnv`.
- `clean` — deletes everything under `--output_dir` except `figures/` and `tables/`.

## External tools

The pipeline shells out to CNVkit twice; both calls are issued from `main.py` with `ulimit -s 65536`:

- **Step 3** — `cnvlib.cnvkit.main()` invoked via an inline `python -c` that monkey-patches `DataFrame.iteritems` / `Series.iteritems` to `.items` for pandas 2.x compatibility, then runs `import-rna -f counts -g <gene_info> -c <corr_file> --output-dir cnr -o cnr/output.txt $(ls *.txt)` from inside `txt/`.
- **Step 5** — `cnvkit.py export cdt *.cnr -o cdt/grpWt.cdt` from inside `wtcnr/`.

## Cross-package data flow

The filesystem under `<output_dir>/` is the contract between stages — each directory has one writer and one or more readers:

| Path | Written by | Read by |
|------|-----------|---------|
| `cluster_exp/Cluster*.parquet`, `cluster_exp/ensembl.csv` | `split_cluster_expr.GeneDataWriter` (Step 1) | `parallel_cluster_runner.ParallelClusterProcessor` and `ClusterMerger` (Step 2); `main.verify_parquet_files` on subsequent runs |
| `txt/<barcode>.txt` | `cluster_processor.{SingleBarcodeProcessor, BarcodesSortedProcessor, SmallClusterProcessor}` (Step 2) | `cnvkit import-rna` subprocess (Step 3) |
| `cluster_summary/<cluster>_barcode_grouping_info.csv` | `cluster_processor.BarcodeMerger` → moved by `OutputDirManager.move_csvs_to_summary` | `call_cnv.ClusterSummary` (Step 6) |
| `cnr/*.cnr` | `cnvkit import-rna` (Step 3) | `weighted_median.WorkflowOrchestrator` → `data_io.CnrReader` (Step 4) |
| `wtcnr/<name>.cnr`, `wtcnr/summary.csv` | `weighted_median.WtMedianWriter` (Step 4) | `cnvkit export cdt` (Step 5) |
| `cdt/grpWt.cdt` | `cnvkit export cdt` (Step 5) | `call_cnv.CNVCallPlotWorkflow` → `data_io.CdtLoader` (Step 6) |
| `tables/*` | `call_cnv.CNVWriter` (Step 6) | end-user / Loupe Browser |
| `figures/*` | `visualization/*` (Steps 1 and 6) | end-user |

Each stage is also re-entrant: `main.py` checks for the presence (and schema, in the case of parquet) of the expected outputs before running a step and skips it if the previous run already produced them.

## Parallelism

Two `multiprocessing.Pool` boundaries exist:

- **Step 2** — one worker per cluster, sized by `--num_processes` (default `cpu_count() - 1`). Each worker owns one cluster's KD-tree merging.
- **Step 4** — one worker per cnr file, sized by `--cores` (default 4). Each worker computes weighted medians for one barcode group's CNR.

Both pools are short-lived and process-isolated; there is no shared mutable state between workers other than the output directory (`wtcnr/summary.csv` is currently written by every Step 4 worker, last-writer-wins — a known issue).

## Dependencies

Runtime dependencies pulled in by `pyproject.toml` and used across the tree: `pandas`, `numpy`, `scipy` (sparse, stats, signal, cluster.hierarchy, spatial.distance), `pyarrow`, `h5py`, `scikit-learn` (`sklearn.neighbors.KDTree`), `wquantiles`, `matplotlib`, `seaborn`, `tqdm`, `cnvkit` / `cnvlib`. Python ≥ 3.9.

## Conventions and design notes

- **Empty `__init__.py`** in every subpackage. There is no re-export layer; importers always use the full module path (`from stmut_hires.weighted_median.parallel_wtmedian import WorkflowOrchestrator`).
- **Static-method classes** are used as namespaces in `visualization/`, `data_io/cdt_loader.py`, `data_io/annotate_reader.py`, and `call_cnv/{save_cnv,cnvscore_permutation,cnvscore_fdr}.py`. The classes never hold state; they exist purely to group related helpers.
- **Configuration is passed, not imported.** Every workflow class takes the same `CNVAnalysisConfig` instance via the constructor. No module reads CLI args or environment variables directly.
- **Reference data is package-bundled.** The three hg38 files in `data/` are resolved with `importlib.resources` so the tool works whether installed via `pip install` or run from a checkout.
- **Re-entrant by design.** Each pipeline stage is guarded by an existence/schema check in `main.py`, so re-running `stmut-hires run` after an interrupted job resumes at the first stage with missing outputs.
