# Changelog

All notable changes to stmut-hires will be documented in this file.

## [2026-01-04]
### Added
- **Sub-command `clean`**: Clean the output folder.
- **Filter Low Quality Barcodes**: Either Automatical predict cutoff (by default) by setting `--bw_method`, or user provided cutoff by setting `--manual_cutoff`

### Fixed
- restructure IO by adding `data_io` module
- fix only pArms genes exist
- automatic add docker tag in CI/CD


## [2026-01-04]
### Added
- **Sub-command `call-cnv`**: Dedicated workflow to rerun CNV calling (Step 6), specifically optimized for X-chromosome gain/loss recovery.
- **Grouping Analytics**: Automated generation of `grouped_spotSummary.csv` and histogram plots to visualize bin grouping statistics.
- **HPC Stability Fix**: Integrated environment activation hooks to resolve `CXXABI_1.3.15` and `libstdc++.so.6` linking issues on enterprise Linux (CentOS 7).
- **Comprehensive Logging**: Replaced standard output prints with a structured logging system for better debugging in headless HPC environments.
- **Unit Testing**: Initial test suite for core functions and CLI entry points.

### Changed
- **Documentation**: Updated README with enterprise-grade hardware specifications and performance benchmarks.


## [2025-12-15]

### Added
- Multiprocessing support for parallel cluster processing
- H5 file format support for gene expression matrices
- Automatic window-size adjustment when gene count < cutoff
- Memory-efficient parquet processing

### Changed
- Improved performance: 11 hours → 3 hours for BD17 dataset
- Updated default cutoff from 500 to 1000 genes
- Enhanced error handling for large datasets

### Fixed
- Memory overflow issues with large expression matrices
- Import path issues in legacy modules
- Window size calculation bug

## [2025-09-21] - Legacy Version

### Added
- Initial implementation of barcode grouping
- CSV-based gene expression processing
- Basic spatial coordinate handling

### Known Issues
- Very slow processing for big cluster data(11+ hours for large datasets)
- Memory issues with datasets > 10GB
- No parallel processing support
- No automate window-size adjustment



