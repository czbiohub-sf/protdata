import warnings
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import pandas as pd


def read_spectronaut(
    file: str,
    intensity_columns: Optional[List[str] | str] = ["PG.Quantity"],
    index_column: Optional[str] = "PG.ProteinGroups",
    sample_column: Optional[str] = "R.FileName",
    sep="\t",
):
    """
    Load Spectronaut results into an AnnData object.

    Parameters
    ----------
    file
        Path to the Spectronaut results file.
    intensity_columns
        Name of the intensity column.
    index_column
        Name of the column to use as protein index.
    sample_column
        Name of the column to use as sample index.
    sep
        File separator.

    Returns
    -------
    :class:`anndata.AnnData` object with:

        - ``X``: intensity matrix (samples x proteins)
        - ``var``: protein metadata (indexed by protein group IDs)
        - ``obs``: sample metadata (indexed by sample names)
    """

    if isinstance(intensity_columns, str):
        intensity_columns = [intensity_columns]
    # Check that intensity columns are all at the same level as the index column
    index_level = index_column.split(".")[0]
    intensity_levels = [col.split(".")[0] for col in intensity_columns]
    if not all(level == index_level for level in intensity_levels):
        raise ValueError(
            f"Intensity columns {intensity_columns} are not all at the same level as the index column {index_column}"
        )

    sample_levels = ("E", "R")
    possible_levels = ("PG", "PEP", "EG", "FG", "F")
    df = _read_csv_auto_decimal(file, sep=sep)

    sample_level = sample_column.split(".")[0]
    sample_idx = sample_levels.index(sample_level)
    var_levels = sample_levels[: sample_idx + 1]
    var_levels = tuple(f"{level}." for level in var_levels)
    var_cols = [
        col
        for col in df.columns
        if col.startswith(var_levels) and df[col].notna().sum() > 0
    ]

    dfp = df.pivot_table(
        index=index_column,
        columns=var_cols,
        values=intensity_columns,
        aggfunc="mean",
        observed=True,
        dropna=True,
    ).T
    dfpx = dfp.loc[intensity_columns[0]]

    # Report if intensities were not unique for pivoting and aggregated. This might not be what the user wants.
    for intensity_column in intensity_columns:
        if df[intensity_column].nunique() > dfpx.size:
            warnings.warn(
                f"{intensity_column} is not unique within {index_column} and {sample_column}, pivoting with mean aggregation. Make sure this is what you want!"
            )

    # Construct obs
    obs = dfpx.index.to_frame().set_index(sample_column)

    # Construct var
    # Only columns that are coarser than the value column make sense to keep in var
    var_levels = tuple(
        f"{level}."
        for level in possible_levels[: possible_levels.index(index_level) + 1]
    )
    varg = df.loc[:, df.columns.str.startswith(var_levels)].groupby(
        index_column, observed=True, dropna=False
    )
    uniquecols = varg.nunique().eq(1).all(axis=0)
    uniquecols = uniquecols[uniquecols].index
    var = varg[uniquecols].first().reindex(dfpx.columns)

    # Get the layers
    layers = {}
    for intensity_column in intensity_columns[1:]:
        layers[intensity_column] = dfp.loc[intensity_column]

    # Create AnnData
    uns = {
        "RawInfo": {
            "Search_Engine": "Spectronaut",
        },
    }
    adata = ad.AnnData(X=dfpx.to_numpy(), obs=obs, var=var, uns=uns, layers=layers)
    return adata


def _read_csv_auto_decimal(
    path: str,
    *,
    sample_rows: int = 2000,
    trials: List[Dict[str, Optional[str]]] = [
        {
            "decimal": ".",
        },
        {
            "decimal": ",",
        },
    ],
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Read a CSV by trying multiple (decimal, thousands) configurations on a sample
    using a fixed `sep` (no separator inference). Picks the best config and
    re-reads the full file once.

    Parameters
    ----------
    path : str
        CSV file path.
    sample_rows : int
        Number of rows to sample for scoring.
    trials : list of dict
        List of {'decimal': <'.' or ','>, 'thousands': <'.' or ',' or None>} to try.
    **kwargs : Any
        Passed through to `pd.read_csv`.

    Returns
    -------
    DataFrame
    """

    scores: List[Tuple[int, Dict[str, Optional[str]]]] = []
    base_kwargs = dict(kwargs)

    for cfg in trials:
        try:
            samp = pd.read_csv(path, nrows=sample_rows, **base_kwargs, **cfg)
            num = samp.select_dtypes(include="number")
            # Score = (non-NaN numeric cells)*10 + (# numeric cols)
            frac_valid = 1.0 - num.isna().to_numpy().mean()
            score = frac_valid * 10000 + num.shape[1]  # or any large scale

            scores.append((score, cfg))
        except Exception:
            scores.append((-1, cfg))

    best_score, best_cfg = max(scores, key=lambda x: x[0])
    if best_score < 0:
        # All trials failed; fall back to single read with caller's sep and defaults
        return pd.read_csv(path, **base_kwargs)

    return pd.read_csv(path, **base_kwargs, **best_cfg)
