import os
import warnings

import numpy as np
import pandas as pd
import pytest

from protdata.io.spectronaut_loader import read_spectronaut

# Absolute path provided by the user for local validation
LOCAL_DATA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../data/spectronaut_normal.tsv")
)


def test_spectronaut_pivot_warning_and_shape():
    with pytest.warns(UserWarning, match="pivoting with mean aggregation"):
        adata = read_spectronaut(
            LOCAL_DATA_PATH,
            index_column="EG.Workflow",
            intensity_columns=["EG.Cscore", "EG.Qvalue"],
        )
    assert adata.shape == (9, 1)


def test_spectronaut_no_warning_and_shape():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        adata = read_spectronaut(
            LOCAL_DATA_PATH,
            index_column="EG.ModifiedSequence",
            intensity_columns=["EG.Cscore", "EG.Qvalue"],
        )
        assert len(caught) == 0
    assert adata.shape == (9, 14)


def _count_nonnull_grouped_values(
    df: pd.DataFrame, index_col: str, sample_col: str, intensity_col: str
) -> int:
    grouped = (
        df[[index_col, sample_col, intensity_col]]
        .groupby([index_col, sample_col], observed=True, dropna=False)[intensity_col]
        .apply(lambda s: pd.notna(s).any())
    )
    return int(grouped.sum())


def _count_nonnull_in_matrix(matrix: np.ndarray) -> int:
    return int(np.count_nonzero(~np.isnan(matrix)))


def test_spectronaut_pivot_nonnull_counts_match_local_data():
    assert os.path.isfile(LOCAL_DATA_PATH)
    # Read raw for ground truth counting (file uses comma decimals)
    df = pd.read_csv(LOCAL_DATA_PATH, sep="\t", decimal=",")
    adata = read_spectronaut(LOCAL_DATA_PATH)

    grouped_nonnull = _count_nonnull_grouped_values(
        df,
        index_col="PG.ProteinGroups",
        sample_col="R.FileName",
        intensity_col="PG.Quantity",
    )
    matrix_nonnull = _count_nonnull_in_matrix(adata.X)

    assert matrix_nonnull == grouped_nonnull


def test_spectronaut_layers_shape_and_presence():
    assert os.path.isfile(LOCAL_DATA_PATH)
    adata = read_spectronaut(
        LOCAL_DATA_PATH,
        intensity_columns=["PG.Quantity", "PG.Cscore"],
    )
    assert "PG.Cscore" in adata.layers
    assert adata.layers["PG.Cscore"].shape == adata.X.shape


def test_spectronaut_uns_and_obs_metadata():
    assert os.path.isfile(LOCAL_DATA_PATH)
    adata = read_spectronaut(LOCAL_DATA_PATH)
    assert adata.uns["RawInfo"]["Search_Engine"] == "Spectronaut"
    # obs index is the sample identifier column
    assert adata.obs.index.name == "R.FileName"
    # expect at least one other sample-level column to be present
    assert "R.Replicate" in adata.obs.columns or "R.Condition" in adata.obs.columns


def test_spectronaut_level_mismatch_raises():
    assert os.path.isfile(LOCAL_DATA_PATH)
    with pytest.raises(ValueError):
        read_spectronaut(
            LOCAL_DATA_PATH,
            index_column="EG.Workflow",
            intensity_columns=["PG.Quantity"],
        )
