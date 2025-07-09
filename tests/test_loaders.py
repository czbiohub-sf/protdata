import os
import pytest
import anndata as ad
from protdata.io.maxquant_loader import read_maxquant
from protdata.io.fragpipe_loader import read_fragpipe
from protdata.io.mztab_loader import read_mztab

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))


@pytest.mark.parametrize(
    "filename,loader",
    [
        ("proteinGroups.txt", read_maxquant),
        ("combined_protein.tsv", read_fragpipe),
        ("SILAC_SQ.mzTab", read_mztab),
    ],
)
def test_loader(filename, loader, tmp_path):
    path = os.path.join(data_dir, filename)
    if not os.path.isfile(path):
        pytest.skip(f"Test data file {filename} not found.")
    adata = loader(path)
    assert isinstance(adata, ad.AnnData)
    assert adata.shape[0] > 0 and adata.shape[1] > 0
    # Make sure the anndata file can be saved and loaded
    test_file = tmp_path / "test.h5ad"
    adata.write_h5ad(test_file)
    adata2 = ad.read_h5ad(test_file)
    assert adata2.shape == adata.shape
    assert (
        adata2.uns["RawInfo"]["Search_Engine"] == adata.uns["RawInfo"]["Search_Engine"]
    )
    if "filter_columns" in adata.uns["RawInfo"]:
        assert list(adata2.uns["RawInfo"]["filter_columns"]) == list(
            adata.uns["RawInfo"]["filter_columns"]
        )
    assert adata2.var.equals(adata.var)
    assert adata2.obs.equals(adata.obs)
