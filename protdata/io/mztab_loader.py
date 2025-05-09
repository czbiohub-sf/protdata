import pandas as pd
import numpy as np
import anndata as ad
from typing import Union
from pyteomics import mztab


def read_mztab(
    file: Union[str, pd.DataFrame],
    index_column: str = "accession",
) -> ad.AnnData:
    """
    Load mzTab protein table into an AnnData object.

    Args:
        file: Path to mzTab file or a pandas DataFrame (protein table).
        intensity_column_prefix: Prefix for intensity columns (default: 'protein_abundance_').
        index_column: Column indicating the protein groups (default: 'accession').

    Returns:
        AnnData object with:
            - X: intensity matrix (proteins x samples)
            - var: protein metadata
            - obs: sample metadata
    """
    if isinstance(file, pd.DataFrame):
        df = file.copy()
    else:
        tables = mztab.MzTab(file)
        df = tables.protein_table

    # Find intensity columns
    intensity_cols = [
        col for col in df.columns if col.startswith("protein_abundance_study_variable[")
    ]
    if not intensity_cols:
        raise ValueError(
            f"No columns starting with 'protein_abundance_study_variable' found."
        )

    # Extract sample names from intensity columns
    sample_names = [col.split("[")[1].split("]")[0] for col in intensity_cols]

    # Build X matrix (proteins x samples)
    X = df[intensity_cols].to_numpy(dtype=np.float32).T

    # Build var (proteins)
    var = df.drop(columns=intensity_cols).copy()
    var.index = df[index_column].astype(str)

    # Build obs (samples)
    obs = pd.DataFrame.from_dict(tables.study_variables, orient="index")
    obs.index.name = "sample"
    obs.index = obs.index.astype(str)

    # Build uns
    uns = {"Search_Engine": df.search_engine.iloc[0]}

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var, uns=uns)
    return adata
