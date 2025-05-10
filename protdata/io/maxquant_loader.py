import pandas as pd
import numpy as np
import anndata as ad
from typing import Optional, Union, List
import warnings


def read_maxquant(
    file: Union[str, pd.DataFrame],
    intensity_column_prefixes: List[str] | str = ["LFQ intensity ", "Intensity "],
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

    Args:
        file: Path to proteinGroups.txt or a pandas DataFrame.
        intensity_column_prefix: Prefix for intensity columns (default: 'LFQ intensity ').
        index_column: Column name for protein IDs (default: 'Protein IDs').
        gene_names_column: Column name for gene names (default: 'Gene names').
        sep: File separator (default: tab).

    Returns:
        AnnData object with:
            - X: intensity matrix (proteins x samples)
            - var: protein metadata
            - obs: sample metadata
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
                layers[prefix.strip()] = df[prefix_cols].to_numpy(dtype=np.float32).T
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

    # Build obs (samples)
    obs = pd.DataFrame(index=sample_names)

    # Build uns
    uns = {"Search_Engine": "MaxQuant"}

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers, uns=uns)
    return adata
