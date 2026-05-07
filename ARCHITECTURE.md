# Software Design Document

## Overview
stmut-hires is a command-line tool for Copy Number Variation (CNV) analysis in high-resolution spatial transcriptomics data. It processes gene expression data from platforms like VisiumHD and StereoSeq, performing clustering, expression aggregation, and CNV calling.

## Architecture

### High-Level Architecture
The system follows a modular pipeline architecture with the following main components:

1. **CLI Interface** (`cli_parser.py`): Parses command-line arguments and subcommands.
2. **Configuration Management** (`config.py`): Handles all configuration with inheritance-based classes.
3. **Main Orchestrator** (`main.py`): Coordinates the entire pipeline execution.
4. **Processing Modules**:
   - `split_cluster_expr/`: Handles cluster expression splitting and processing.
   - `parallel/`: Manages parallel processing and output management.
   - `weighted_median/`: Performs weighted median calculations.
   - `call_cnv/`: Executes CNV calling and visualization.
5. **Data Package** (`data/`): Contains reference data files.

### Pipeline Flow
The full pipeline consists of 6 steps:

1. **Cluster Expression Processing**: Splits gene expression data by clusters and saves to Parquet files.
2. **Parallel Barcode Merging**: Processes barcodes in parallel, merging spatial and expression data.
3. **CNR File Generation**: Uses CNVkit to generate Copy Number Ratio files.
4. **Weighted Median Calculation**: Applies weighted median smoothing to CNR data.
5. **CDT Export**: Exports data to CDT format for visualization.
6. **CNV Calling and Plotting**: Performs final CNV analysis and generates plots.

### Key Design Patterns

#### Configuration Inheritance
- `BaseConfig`: Abstract base with common settings.
- `InitialStepConfig`: Adds input file paths.
- `MergerConfig`: Adds spatial file handling.
- `CNVAnalysisConfig`: Adds CNV-specific parameters.

#### Modular Processing
Each major step is encapsulated in its own module with clear responsibilities:
- `ClusterExpressionProcessor`: Handles step 1.
- `ParallelClusterProcessor`: Handles step 2.
- `WorkflowOrchestrator`: Handles step 4.
- `CNVCallPlotWorkflow`: Handles step 6.

#### Parallel Processing
- Uses multiprocessing for CPU-intensive tasks.
- Configurable number of processes and cores.
- Memory-efficient processing of large datasets.

### Data Flow
```
Input Files (H5, CSV, Parquet)
    ↓
Cluster Processing (Step 1)
    ↓
Parallel Merging (Step 2)
    ↓
CNR Generation (Step 3 - External: CNVkit)
    ↓
Weighted Median (Step 4)
    ↓
CDT Export (Step 5 - External: CNVkit)
    ↓
CNV Calling & Visualization (Step 6)
    ↓
Output Files (Plots, CSVs, etc.)
```

### Dependencies
- **Core**: pandas, pyarrow, numpy, scipy
- **Visualization**: seaborn, matplotlib
- **External Tools**: CNVkit (via subprocess)
- **Parallel**: multiprocessing, concurrent.futures

### Error Handling
- File existence validation in config classes.
- Subprocess error checking for external tools.
- Logging throughout the pipeline.

### Performance Considerations
- Parquet format for intermediate data storage.
- Parallel processing for scalability.
- Memory-efficient data structures.
- Benchmarking shows linear scaling with data size.

## Future Improvements
- Add unit test coverage.
- Implement configuration file support (YAML).
- Add Docker containerization improvements.
- Enhance error handling and recovery.
- Add type hints throughout codebase.
- Implement CI/CD pipeline.</content>
<parameter name="filePath">/analysis_directory/limin/apps/stmut-hires/ARCHITECTURE.md