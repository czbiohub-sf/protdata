# protdata

[![Test](https://github.com/czbiohub-sf/protdata/actions/workflows/test.yml/badge.svg)](https://github.com/czbiohub-sf/protdata/actions/workflows/test.yml)
![PyPI - Version](https://img.shields.io/pypi/v/protdata)
[![docs online](https://img.shields.io/badge/docs-online-blue)](https://protdata-czbiohub.vercel.app/)

Proteomics data loaders for the [AnnData](https://anndata.readthedocs.io/) format.

This package provides loader functions to import proteomics data (e.g., MaxQuant) into the AnnData structure for downstream analysis and integration with single-cell and multi-omics workflows.

## Features

- **Multiple formats**: Support for MaxQuant, FragPipe, DIA-NN, and mzTab files
- **Reads metadata**: Automatically extracts and organizes sample and protein metadata

## Installation

Protdata has minimal dependencies listed in [pyproject.toml](https://github.com/czbiohub-sf/protdata/blob/main/pyproject.toml)

To install the latest release from PyPI, run:

```bash
pip install protdata
```

Or install from source:
```bash
git clone https://github.com/czbiohub-sf/protdata.git
cd protdata
pip install -e . # or make setup-develop for developers
```

## Usage Example

### MaxQuant Import

You can download an example proteinGroups [file here](https://github.com/czbiohub-sf/protdata/raw/main/data/proteinGroups.txt)
```python
import protdata

adata = load_maxquant_to_anndata("/path/to/proteinGroups.txt")
print(adata)
```

### DIA-NN Import

You can download an example DIA-NN report [file here](https://github.com/czbiohub-sf/protdata/raw/main/data/report.pg_matrix.tsv)

```python
from protdata.io import read_diann

adata = read_diann("/path/to/report.pg_matrix.tsv")
print(adata)
```

### FragPipe Import

You can download an example FragPipe output [file here](https://github.com/czbiohub-sf/protdata/raw/main/data/combined_protein.tsv)

```python
from protdata.io import read_fragpipe

adata = read_fragpipe("/path/to/combined_protein.tsv")
print(adata)
```

### mzTab Import

You can download an example mzTab [file here](https://github.com/czbiohub-sf/protdata/raw/main/data/SILAC_SQ.mzTab)

```python
from protdata.io import read_mztab

adata = read_mztab("/path/to/SILAC_SQ.mzTab")
print(adata)
```
## Authors

`protdata` is created and maintained by the [Computational Biology Platform](https://www.czbiohub.org/comp-biology/) at the [Chan Zuckerberg Biohub San Francisco](https://www.czbiohub.org/sf/).

To get in touch please use the [GihHub issues](https://github.com/czbiohub-sf/protdata/issues) page.

## Contributing

If you want to contribute to `protdata`, please read the [Contribution Guide](https://prodata.sf.czbiohub.org/contributing)

## Changelog
See [Release Notes](https://prodata.sf.czbiohub.org/release_notes)
