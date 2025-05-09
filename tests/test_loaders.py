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
def test_loader(filename, loader):
    path = os.path.join(data_dir, filename)
    if not os.path.isfile(path):
        pytest.skip(f"Test data file {filename} not found.")
    adata = loader(path)
    assert isinstance(adata, ad.AnnData)
    assert adata.shape[0] > 0 and adata.shape[1] > 0
