import warnings
from pathlib import Path

import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import biotapy as bt
from biotapy.datasets import _biocrust, _remote

# Synthetic tables in the layout of mmvec's soil example: the metabolite table lacks one microbe sample and lists the
# others in another order. No data of the example is copied.
_MICROBES = AnnData(
    X=sp.csr_matrix(np.array([[10, 0, 3], [4, 6, 0], [0, 2, 8]])),
    obs=pd.DataFrame(index=["3min_early", "3min_late", "9hr_early"]),
    var=pd.DataFrame(index=["rplo 1 (Cyanobacteria)", "rplo 2 (Firmicutes)", "rplo 3 (Proteobacteria)"]),
)
_METABOLITES = AnnData(
    X=sp.csr_matrix(np.array([[1500.5, 2.0], [800.25, 31.75]])),
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
    np.testing.assert_array_equal(mdata["metabolites"].X.toarray(), [[800.25, 31.75], [1500.5, 2.0]])


def test_fetches_both_pinned_files_before_reading(fetched):
    bt.datasets.biocrust()
    assert fetched == ["biocrust_microbes.biom", "biocrust_metabolites.biom"]


@pytest.mark.parametrize("page", ["guide/datasets.md", "tutorials/multiomics.md"])
def test_the_pages_link_mmvec_s_example_at_the_pinned_commit(page):
    commit = _remote._MMVEC.split("/mmvec/")[1].split("/")[0]
    text = (Path(__file__).parents[2] / "docs" / page).read_text(encoding="utf-8")
    assert f"https://github.com/biocore/mmvec/tree/{commit}/examples/soils" in text


@pytest.mark.network
def test_biocrust_downloads_and_loads():
    mdata = bt.datasets.biocrust()
    assert {key: mod.shape for key, mod in mdata.mod.items()} == {"taxa": (19, 466), "metabolites": (19, 85)}
    assert "9hr_late" not in mdata.obs_names


@pytest.mark.network
def test_biocrust_prose_figures_hold():
    # The figures the guide, the tutorial and the docstring quote about the downloaded files.
    paths = [Path(_biocrust._fetch(name)) for name in _biocrust.FILES]
    assert round(sum(path.stat().st_size for path in paths) / 1000) == 135
    microbes = bt.io.read_biom(paths[0])
    assert microbes.n_obs == 20
    assert microbes.var_names[np.asarray(microbes.X.sum(axis=0)).ravel().argmax()] == "rplo 1 (Cyanobacteria)"
