Release Notes
=============

This page contains information about protdata releases and changes.

Version 0.1.5
--------------

**Features**

- Added support for Spectronaut reports with :func:`protdata.io.read_spectronaut`
- Added Spectronaut example data and tests

**Bug fixes**

- Fixed h5ad serialization of MaxQuant data: forward slashes (``/``) in column names
  are replaced with underscores, since hdf5 interprets them as groups
- Fixed links in the README

**Breaking changes**

- :func:`protdata.io.read_maxquant` now replaces ``/`` with ``_`` in all column names,
  so ``var`` columns that contained a ``/`` are renamed

Version 0.1.4
--------------

**Features**

- Added documentation to the project
- Added tests
- Updated links to example files in the README

Version 0.1.3
--------------

**Features**

- Added support for mzTab files with :func:`protdata.io.read_mztab`
- Improved MaxQuant loader with better layer support
- Added FragPipe combined_protein.tsv support
- DIA-NN report.pg_matrix.tsv support

**Bug fixes**

- Fixed h5ad serialization issues with mixed data types
- Improved error handling for missing columns
- Better string conversion for metadata

**Breaking changes**

- None

Version 0.1.2
--------------

**Features**

- Initial MaxQuant proteinGroups.txt support
- Basic AnnData structure implementation

**Bug fixes**

- None

**Breaking changes**

- None

Version 0.1.1
--------------

**Features**

- Initial release
- Basic project structure

**Bug fixes**

- None

**Breaking changes**

- None 