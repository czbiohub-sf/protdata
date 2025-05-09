import pandas as pd
import numpy as np
import anndata as ad
from typing import Union
import os


def read_diann(
    file: Union[str, pd.DataFrame],
    index_column: str = "Protein.Group",
    sep: str = "\t",
) -> ad.AnnData:
    """
    Load DIA-NN protein group matrix (report.pg_matrix.tsv) into an AnnData object.

    Args:
        file: Path to DIA-NN report.pg_matrix.tsv or a pandas DataFrame.
        intensity_suffix: Suffix for intensity columns (default: '_Intensity').
        index_column: Column name for protein group IDs (default: 'Protein.Group').
        sep: File separator (default: tab).

    Returns:
        AnnData object with:
            - X: intensity matrix (proteins x samples)
            - var: protein metadata
            - obs: sample metadata
    """
    if isinstance(file, pd.DataFrame):
        df = file.copy()
    else:
        df = pd.read_csv(file, sep=sep, low_memory=False)

    # Columns that are not sample columns (from DIANNLoader.no_sample_column)
    no_sample_column = [
        "PG.Q.value",
        "Global.PG.Q.value",
        "PTM.Q.value",
        "PTM.Site.Confidence",
        "PG.Quantity",
        "Protein.Group",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "First.Protein.Description",
        "contamination_library",
    ]

    df.columns = [
        os.path.basename(col)
        if isinstance(col, str) and col not in no_sample_column
        else str(col)
        for col in df.columns
    ]

    # Add contamination column if not present
    if "contamination_library" not in df.columns:
        df["contamination_library"] = False

    # Find intensity columns (all columns not in no_sample_column)
    intensity_cols = df.columns.difference(no_sample_column, sort=False).tolist()
    if not intensity_cols:
        raise ValueError(f"No intensity columns found.")

    # Build X matrix (proteins x samples)
    X = df[intensity_cols].to_numpy(dtype=np.float32).T

    # Build var (proteins)
    var = df[df.columns.intersection(no_sample_column, sort=False)].copy()
    var.index = df[index_column]

    # Build obs (samples)
    obs = pd.DataFrame(index=intensity_cols)

    # Build uns
    uns = {"Search_Engine": "DIANN"}

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns = uns
    return adata
