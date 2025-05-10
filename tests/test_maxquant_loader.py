import pandas as pd
import numpy as np
from protdata.io.maxquant_loader import read_maxquant


def test_load_maxquant_to_anndata():
    data = {
        "Protein IDs": ["P1", "P2"],
        "Gene names": ["G1", "G2"],
        "LFQ intensity SampleA": [100.0, 200.0],
        "LFQ intensity SampleB": [150.0, 250.0],
    }
    df = pd.DataFrame(data)
    adata = read_maxquant(df, intensity_column_prefixes=["LFQ intensity "])
    assert adata.shape == (2, 2)  # 2 proteins x 2 samples
    np.testing.assert_array_equal(
        adata.X, np.array([[100.0, 150.0], [200.0, 250.0]], dtype=np.float32).T
    )
    assert list(adata.var.index) == ["P1", "P2"]
    assert list(adata.obs.index) == ["SampleA", "SampleB"]
    assert "Gene names" in adata.var.columns
