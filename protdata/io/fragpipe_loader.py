import pandas as pd
import numpy as np
import anndata as ad
from typing import Union
import warnings


def read_fragpipe(
    file: Union[str, pd.DataFrame],
    intensity_column_suffixes=[
        " MaxLFQ Intensity",
        " Spectral Count",
        " Unique Spectral Count",
    ],
    index_column: str = "Protein ID",
    sep: str = "\t",
) -> ad.AnnData:
    """
    Load FragPipe-Philosopher protein group matrix into an AnnData object.

    Args:
        file: Path to combined_protein.tsv or a pandas DataFrame.
        intensity_column_prefix: Prefix for intensity columns (default: '[sample] MaxLFQ Intensity ').
        index_column: Column name for protein IDs (default: 'Protein').
        gene_names_column: Column name for gene names (default: 'Gene Names').
        confidence_column: Column name for protein probability (default: 'Protein Probability').
        sep: File separator (default: tab).

    Returns:
        AnnData object with:
            - X: intensity matrix (proteins x samples)
            - var: protein metadata
            - obs: sample metadata
    """
    if isinstance(intensity_column_suffixes, str):
        intensity_column_suffixes = [intensity_column_suffixes]
    main_intensity_suffix = intensity_column_suffixes[0]

    if isinstance(file, pd.DataFrame):
        df = file.copy()
    else:
        df = pd.read_csv(file, sep=sep, low_memory=False)

    # Find intensity columns
    intensity_cols = [col for col in df.columns if col.endswith(main_intensity_suffix)]
    if not intensity_cols:
        raise ValueError(f"No columns ending with '{main_intensity_suffix}' found.")

    # Extract sample names from intensity columns
    sample_names = [col[: -len(main_intensity_suffix)] for col in intensity_cols]

    # Build X matrix (proteins x samples)
    X = df[intensity_cols].to_numpy(dtype=np.float32).T
    # If there are more intensity suffixes we store them as layers
    layers = {}
    if len(intensity_column_suffixes) > 1:
        for suffix in intensity_column_suffixes[1:]:
            suffix_cols = [
                col
                for col in df.columns
                if col in [sample_name + suffix for sample_name in sample_names]
            ]
            if len(suffix_cols) == len(sample_names):
                layers[suffix.strip()] = df[suffix_cols].to_numpy(dtype=np.float32).T
            else:
                warnings.warn(
                    f"Number of columns for suffix '{suffix}' does not match number of samples."
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
    uns = {"Search_Engine": "FragPipe_Philosopher"}

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers)
    adata.uns = uns

    return adata
