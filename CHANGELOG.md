# Changelog

All notable changes to stmut-hires will be documented in this file.

## [Unreleased]

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



