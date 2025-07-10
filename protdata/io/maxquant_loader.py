from __future__ import annotations

import warnings
from typing import List, Union

import anndata as ad
import numpy as np
import pandas as pd

from .utils import cleanup_obsvar


def read_maxquant(
    file: Union[str, pd.DataFrame],
    intensity_column_prefixes: List[str] | str = [
        "LFQ intensity ",
        "Intensity ",
        "MS/MS count ",
    ],
    index_column: str = "Protein IDs",
    filter_columns: list[str] = [
        "Only identified by site",
        "Reverse",
        "Potential contaminant",
    ],
    # gene_names_column: str = "Gene names",
    sep: str = "\t",
) -> ad.AnnData:
    """
    Load MaxQuant proteinGroups.txt into an AnnData object.

    Parameters
    ----------
    file
        Path to the MaxQuant proteinGroups.txt file or a pandas DataFrame containing the data.
    intensity_column_prefixes
        Prefix(es) for intensity columns to extract.
        The first prefix is used for the main matrix (X), others are stored as layers if present.
    index_column
        Column name to use as protein index.
    filter_columns
        Columns to use for filtering out contaminants or unwanted entries.
    sep
        File separator if reading from file.

    Returns
    -------
    anndata.AnnData
        AnnData object with:

        - X: intensity matrix (samples x proteins)
        - var: protein metadata (indexed by protein IDs)
        - obs: sample metadata (indexed by sample names)
        - layers: additional intensity matrices if multiple intensity column prefixes are provided
    """
    if isinstance(intensity_column_prefixes, str):
        intensity_column_prefixes = [intensity_column_prefixes]

    main_intensity_column = intensity_column_prefixes[0]
    if isinstance(file, pd.DataFrame):
        df = file
    else:
        df = pd.read_csv(file, sep=sep, low_memory=False)

    # Find intensity columns
    intensity_cols = [
        col for col in df.columns if col.startswith(main_intensity_column)
    ]
    if not intensity_cols:
        raise ValueError(f"No columns starting with '{main_intensity_column}' found.")

    # Extract sample names from intensity columns
    sample_names = [col[len(main_intensity_column) :] for col in intensity_cols]

    # Build X matrix (proteins x samples)
    X = df[intensity_cols].to_numpy(dtype=np.float32).T

    # If there are more intensity suffixes we store them as layers
    layers = {}
    if len(intensity_column_prefixes) > 1:
        for prefix in intensity_column_prefixes[1:]:
            prefix_cols = [
                col
                for col in df.columns
                if col in [prefix + sample_name for sample_name in sample_names]
            ]
            if len(prefix_cols) == len(sample_names):
                # Cannot have '/' in the key (hdf5 interprets as a group for serialization)
                layers[prefix.strip().replace("/", "_")] = (
                    df[prefix_cols].to_numpy(dtype=np.float32).T
                )
            else:
                warnings.warn(
                    f"Number of columns for prefix '{prefix}' does not match number of samples."
                )

    # Build var (proteins)
    # A metadata column is anything that does not contain a sample name
    sample_columns = np.array(
        [any(sample_name in col for sample_name in sample_names) for col in df.columns]
    )
    var = df.loc[:, ~sample_columns].copy()
    var.index = df[index_column]
    var = cleanup_obsvar(var)

    # Build obs (samples)
    obs = pd.DataFrame(index=sample_names)
    obs = cleanup_obsvar(obs)

    # Build uns
    uns = {
        "RawInfo": {
            "Search_Engine": "MaxQuant",
            "filter_columns": filter_columns,
        },
    }

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers, uns=uns)
    return adata
