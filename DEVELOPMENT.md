# Development Status

## Current Work (Sep 2025)

### Finished.
- [ ] **Sub-command to rerun X-chr gainLoss**
  - Statue: 100% complete
  - Issue: it is typically difficult to infer copy number alterations on the X-chr with gene expression data. One copy of the X-chr is silenced via X-inactivation. If the inactive copy of X is subjected to a CNA, it would not show up in the gene expression data. If the active copy is gained, better to include X-chr just to rerun the RankedBySimilarity analysis.
  - Next: add a sub-command to only run step3  by taking X-chr gain or loss.
  - Files: 

- [ ] **Plot grouping statistic as standard output**
  - Status: 100% complete
  - Issue: Missing statistic of how many bins to form a new bin
  - Next: Generate a histogram to visualize grouping statistics; grouped_spotSummary.csv (.1 not -1 in barcode); output --caseCNV_summ.csv
  - Files


### Planned
- [ ] **Take Stereo-seq data**
  - Priority: High
  - Description: add a new class to infer CNV from Stereo-seq.
  - `local`-level merge statistic file/histogram plot.
  - Timeline: Dec 2026

- [ ] **Configuration file support**
  - Priority: Low
  - Description: YAML config files instead of command-line args
  - Timeline: Q1 2025


## Performance Benchmarks

### BD17 Dataset (4M cells)
| Version | Processing Time | Memory Usage | Notes |
|---------|----------------|--------------|-------|
| Legacy | 11 hours | 32GB | CSV processing |
| Current | 3 hours | 16GB | Parquet + multiprocessing |
| Target | 1 hour | 8GB | GPU acceleration |

### Test Datasets
- **Small**: 100K cells - Legacy: 30min, Current: 5min
- **Medium**: 1M cells - Legacy: 3 hours, Current: 45min
- **Large**: 4M cells - Legacy: 11 hours, Current: 3 hours


## Code Quality

### Current Issues
- [ ] Add unit tests for core functions
- [ ] Improve error handling in multiprocessing
- [ ] Add logging instead of print statements
- [ ] Document all function parameters

### Code Coverage
- Core functions: 60%
- Module functions: 40%
- CLI functions: 20%

## Next Release Goals

### v2.1 (Target: December 15, 2024)
- [ ] Fix auto window-size adjustment
- [ ] Add comprehensive error handling
- [ ] Improve documentation
- [ ] Add example datasets

### v2.2 (Target: January 15, 2025)
- [ ] GPU acceleration support
- [ ] Configuration file support
- [ ] Performance optimizations
- [ ] Docker container

## Development Notes

### Git Workflow
- `main` branch: Stable releases
- `develop` branch: Integration branch
- `feature/*` branches: New features
- `hotfix/*` branches: Bug fixes

### Testing Strategy
- Manual testing with BD17 dataset
- Performance regression testing
- Memory usage monitoring
- Cross-platform compatibility (Linux/Mac)

### Dependencies
- Python 3.9+
- pandas, numpy, scikit-learn
- dask, pyarrow
- h5py (for H5 support)
- multiprocessing (built-in)

## Known Issues

### High Priority
- Window size auto-adjustment not working
- Memory spikes during large dataset processing
- No progress bar for long-running operations

### Medium Priority
- Import path issues in some environments
- Limited error messages for debugging
- No validation of input file formats

### Low Priority
- Documentation could be more detailed
- No automated testing
- Hard-coded parameters in some functions

