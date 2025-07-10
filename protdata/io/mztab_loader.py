from typing import Union

import anndata as ad
import numpy as np
import pandas as pd
from pyteomics import mztab

from .utils import cleanup_obsvar


def read_mztab(
    file: Union[str, pd.DataFrame],
    index_column: str = "accession",
) -> ad.AnnData:
    """
    Load mzTab protein table into an AnnData object.

    Parameters
    ----------
    file
        Path to mzTab file or a pandas DataFrame containing the protein table.
    index_column
        Column indicating the protein groups.

    Returns
    -------
    anndata.AnnData
        AnnData object with:

        - X: intensity matrix (proteins x samples)
        - var: protein metadata (indexed by protein accession)
        - obs: sample metadata (indexed by sample names)
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
            "No columns starting with 'protein_abundance_study_variable' found."
        )

    # Build X matrix (proteins x samples)
    X = df[intensity_cols].to_numpy(dtype=np.float32).T

    # Build var (proteins)
    var = df.drop(columns=intensity_cols).copy()
    var.index = df[index_column].astype(str)
    var = cleanup_obsvar(var)

    # Build obs (samples)
    obs = pd.DataFrame.from_dict(tables.study_variables, orient="index")
    obs.index.name = "sample"
    obs.index = obs.index.astype(str)
    obs = cleanup_obsvar(obs)

    # Build uns
    uns = {"RawInfo": {"Search_Engine": str(df.search_engine.iloc[0])}}

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var, uns=uns)
    return adata
