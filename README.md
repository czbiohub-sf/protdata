# protdata

[![Test](https://github.com/czbiohub-sf/protdata/actions/workflows/test.yml/badge.svg)](https://github.com/czbiohub-sf/protdata/actions/workflows/test.yml)

Proteomics data loaders for the [AnnData](https://anndata.readthedocs.io/) format.

This package provides loader functions to import proteomics data (e.g., MaxQuant) into the AnnData structure for downstream analysis and integration with single-cell and multi-omics workflows.

## Features
- Load MaxQuant output into AnnData
- Designed for extensibility (other loaders coming soon)

## Installation
```bash
pip install protdata
```

## Usage Example

### MaxQuant Import

You can download an example proteinGroups [file here](https://zenodo.org/records/3774452/files/MaxQuant_Protein_Groups.tabular?download=1)
```python
import protdata

adata = load_maxquant_to_anndata("/path/to/proteinGroups.txt")
print(adata)
``` 

### DIA-NN Import

You can download an example DIA-NN report [file here](https://github.com/vdemichev/DiaNN/raw/master/diann-output-examples/report.pg_matrix.tsv)

```python
from protdata.io import read_diann

adata = read_diann("/path/to/report.pg_matrix.tsv")
print(adata)
```

### FragPipe Import

You can download an example FragPipe output [file here](https://github.com/Nesvilab/philosopher/blob/master/testdata/combined_protein.tsv?raw=true)

```python
from protdata.io import read_fragpipe

adata = read_fragpipe("/path/to/combined_protein.tsv")
print(adata)
```

### mzTab Import

You can download an example mzTab [file here](https://raw.githubusercontent.com/HUPO-PSI/mzTab/refs/heads/master/examples/1_0-Proteomics-Release/SILAC_SQ.mzTab)

```python
from protdata.io import read_mztab

adata = read_mztab("/path/to/SILAC_SQ.mzTab")
print(adata)
```



