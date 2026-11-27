import os

import anndata as ad
import numpy as np
from pyteomics import mztab

from protdata.io.mztab_loader import read_mztab

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))


def test_read_mztab_dataframe_matches_file(tmp_path):
    path = os.path.join(data_dir, "SILAC_SQ.mzTab")
    from_file = read_mztab(path)
    from_df = read_mztab(mztab.MzTab(path).protein_table)

    assert list(from_df.obs_names) == list(from_file.obs_names) == ["1", "2"]
    assert from_df.obs.index.name == "sample"
    assert list(from_df.var_names) == list(from_file.var_names)
    assert from_df.var.equals(from_file.var)
    np.testing.assert_array_equal(from_df.X, from_file.X)
    assert from_file.obs["description"].tolist() == [
        "heat shock response of control",
        "heat shock response of treatment",
    ]

    h5ad_path = tmp_path / "mztab_df.h5ad"
    from_df.write_h5ad(h5ad_path)
    reloaded = ad.read_h5ad(h5ad_path)
    assert reloaded.shape == from_df.shape
    assert list(reloaded.obs_names) == ["1", "2"]
