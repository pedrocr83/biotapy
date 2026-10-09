from contextlib import nullcontext

import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from hypothesis import given
from hypothesis import strategies as st
from mudata import MuData

import biotapy as bt


def _metabolites(samples: list[str]) -> AnnData:
    values = np.arange(len(samples) * 3, dtype=np.float64).reshape(len(samples), 3)
    return AnnData(
        X=sp.csr_matrix(values),
        obs=pd.DataFrame({"batch": ["a"] * len(samples)}, index=samples),
        var=pd.DataFrame(index=["m1", "m2", "m3"]),
    )


def test_aligned_modalities_become_one_mudata():
    tdata = bt.datasets.toy()
    mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata.obs_names.tolist())})
    assert isinstance(mdata, MuData)
    assert list(mdata.mod) == ["taxa", "metabolites"]
    assert mdata.obs_names.tolist() == tdata.obs_names.tolist()
    assert mdata["metabolites"].obs_names.equals(mdata["taxa"].obs_names)


def test_keeps_the_samples_every_modality_has_in_the_first_ones_order():
    tdata = bt.datasets.toy()
    with pytest.warns(UserWarning, match=r"dropped: taxa 3 of 6, metabolites 1 of 4"):
        mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(["s4", "s3", "s2", "s9"])})
    assert mdata.obs_names.tolist() == ["s2", "s3", "s4"]
    for name in ("taxa", "metabolites"):
        assert mdata[name].obs_names.tolist() == ["s2", "s3", "s4"]
    np.testing.assert_array_equal(mdata["metabolites"].X.toarray()[:, 0], [6.0, 3.0, 0.0])


def test_keeps_a_treedata_modality_and_its_tree():
    tdata = bt.datasets.toy()
    mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata.obs_names.tolist())})
    taxa = mdata["taxa"]
    assert type(taxa) is bt._core.TreeData
    assert set(bt._core.tree_tips(bt._core.get_tree(taxa))) == set(taxa.var_names)


def test_function_table_modalities_go_side_by_side():
    table = bt.datasets.toy_humann()
    mdata = bt.io.to_mudata({**table.mod, "taxa": bt.datasets.toy()})
    assert list(mdata.mod) == ["function", "function_by_taxon", "taxa"]


def test_returns_copies_and_keeps_the_input(assert_unchanged):
    tdata = bt.datasets.toy()
    metabolites = _metabolites(["s1", "s2"])
    before = tdata.copy()
    with pytest.warns(UserWarning):
        mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
    assert_unchanged(before, tdata)
    assert not mdata["taxa"].is_view and mdata["taxa"] is not tdata
    mdata["metabolites"].X[0, 1] = 99.0
    assert metabolites.X[0, 1] == 1.0


def test_h5mu_drops_a_treedata_modality_tree(tmp_path):
    # The documented reason to save a tree-bearing modality with write_h5td too (contracts/tree-access).
    mdata = bt.io.to_mudata({"taxa": bt.datasets.toy()})
    mdata.write_h5mu(tmp_path / "study.h5mu")
    assert type(mudata.read_h5mu(tmp_path / "study.h5mu")["taxa"]) is AnnData


def test_single_modality():
    mdata = bt.io.to_mudata({"taxa": bt.datasets.toy()})
    assert list(mdata.mod) == ["taxa"] and mdata.n_obs == 6


def test_no_shared_sample_raises():
    with pytest.raises(ValueError, match="modalities share no sample"):
        bt.io.to_mudata({"taxa": bt.datasets.toy(), "metabolites": _metabolites(["x1", "x2"])})


def test_empty_mapping_raises():
    with pytest.raises(ValueError, match="modalities is empty"):
        bt.io.to_mudata({})


def test_a_mudata_entry_raises_naming_the_fix():
    with pytest.raises(TypeError, match=r"modalities\['function'\] must be an AnnData, got MuData"):
        bt.io.to_mudata({"function": bt.datasets.toy_humann()})


def test_repeated_sample_ids_raise():
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        metabolites = _metabolites(["s1", "s1"])
    with pytest.raises(ValueError, match=r"modalities\['metabolites'\] repeats sample ids: \['s1'\]"):
        bt.io.to_mudata({"taxa": bt.datasets.toy(), "metabolites": metabolites})


@pytest.mark.parametrize("name", ["", 3])
def test_bad_modality_name_raises(name):
    with pytest.raises(TypeError, match="modality names must be non-empty strings"):
        bt.io.to_mudata({name: bt.datasets.toy()})


@given(st.lists(st.sampled_from([f"s{i}" for i in range(1, 9)]), min_size=1, max_size=8, unique=True))
def test_every_modality_holds_the_same_samples(samples):
    tdata = bt.datasets.toy()
    shared = [name for name in tdata.obs_names if name in samples]
    metabolites = _metabolites(samples)
    if not shared:
        with pytest.raises(ValueError, match="share no sample"):
            bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
        return
    with pytest.warns(UserWarning) if len(shared) < max(6, len(samples)) else nullcontext():
        mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
    assert mdata.obs_names.tolist() == shared
    assert mdata["metabolites"].obs_names.tolist() == shared
