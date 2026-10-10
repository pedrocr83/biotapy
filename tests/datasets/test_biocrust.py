import warnings

import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import biotapy as bt
from biotapy.datasets import _biocrust

# Synthetic tables in the layout of mmvec's soil example: the metabolite table lacks one microbe sample and lists the
# others in another order. No data of the example is copied.
_MICROBES = AnnData(
    X=sp.csr_matrix(np.array([[10, 0, 3], [4, 6, 0], [0, 2, 8]])),
    obs=pd.DataFrame(index=["3min_early", "3min_late", "9hr_early"]),
    var=pd.DataFrame(index=["rplo 1 (Cyanobacteria)", "rplo 2 (Firmicutes)", "rplo 3 (Proteobacteria)"]),
)
_METABOLITES = AnnData(
    X=sp.csr_matrix(np.array([[2782242.25, 1.0], [3151923.75, 751234.5]])),
    obs=pd.DataFrame(index=["9hr_early", "3min_early"]),
    var=pd.DataFrame(index=["adenine", "uracil"]),
)


@pytest.fixture
def fetched(tmp_path, monkeypatch):
    # As in test_hmp2.py: the only offline route to the loader is its private _fetch (R11.4).
    bt.io.write_biom(_MICROBES, tmp_path / "biocrust_microbes.biom")
    bt.io.write_biom(_METABOLITES, tmp_path / "biocrust_metabolites.biom")
    names = []
    monkeypatch.setattr(_biocrust, "_fetch", lambda name: names.append(name) or str(tmp_path / name))
    return names


def test_keeps_the_samples_both_tables_have_in_the_microbe_order_without_a_warning(fetched):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mdata = bt.datasets.biocrust()
    assert isinstance(mdata, mudata.MuData) and list(mdata.mod) == ["taxa", "metabolites"]
    for mod in mdata.mod.values():
        assert mod.obs_names.tolist() == ["3min_early", "9hr_early"]


def test_values_and_units_come_from_the_reader(fetched):
    mdata = bt.datasets.biocrust()
    assert mdata["taxa"].uns["biotapy"]["x_kind"] == "counts"
    assert mdata["metabolites"].uns["biotapy"]["x_kind"] == "abundance"
    np.testing.assert_array_equal(mdata["taxa"].X.toarray(), [[10, 0, 3], [0, 2, 8]])
    np.testing.assert_array_equal(mdata["metabolites"].X.toarray(), [[3151923.75, 751234.5], [2782242.25, 1.0]])


def test_fetches_both_pinned_files_before_reading(fetched):
    bt.datasets.biocrust()
    assert fetched == ["biocrust_microbes.biom", "biocrust_metabolites.biom"]


@pytest.mark.network
def test_biocrust_downloads_and_loads():
    mdata = bt.datasets.biocrust()
    assert {key: mod.shape for key, mod in mdata.mod.items()} == {"taxa": (19, 466), "metabolites": (19, 85)}
    assert "9hr_late" not in mdata.obs_names
