# protdata

Proteomics data loaders for the [AnnData](https://anndata.readthedocs.io/) format.

This package provides loader functions to import proteomics data (e.g., MaxQuant) into the AnnData structure for downstream analysis and integration with single-cell and multi-omics workflows.

## Features
- Load MaxQuant output into AnnData
- Designed for extensibility (other loaders coming soon)

## Installation
```bash
pip install .
```

## Usage Example
```python
from protdata.io.maxquant_loader import load_maxquant_to_anndata
adata = load_maxquant_to_anndata("/path/to/proteinGroups.txt")
print(adata)
``` 