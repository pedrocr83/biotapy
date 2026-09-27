import json
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt
from biotapy._core import get_tree


def _adata(dense) -> ad.AnnData:
    return ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )


def _depths(adata) -> list[int]:
    return np.asarray(adata.X.sum(axis=1)).ravel().astype(int).tolist()


# toy() sample depths: s1..s6 = 68, 66, 77, 76, 82, 73.


def test_default_depth_is_the_smallest_sample():
    out = bt.pp.rarefy(bt.datasets.toy(), seed=0)
    assert out.n_obs == 6 and _depths(out) == [66] * 6


def test_sample_at_exactly_depth_is_kept():
    with pytest.warns(UserWarning, match=r"dropped 2 sample\(s\).*\['s1', 's2'\]"):
        out = bt.pp.rarefy(bt.datasets.toy(), depth=73, seed=0)
    assert list(out.obs_names) == ["s3", "s4", "s5", "s6"] and _depths(out) == [73] * 4


def test_no_count_exceeds_the_original():
    tdata = bt.datasets.toy()
    out = bt.pp.rarefy(tdata, depth=50, seed=1)
    assert (out.X.toarray() <= tdata[:, out.var_names].X.toarray()).all()


def test_same_seed_same_result_and_generator_accepted():
    first = bt.pp.rarefy(bt.datasets.toy(), depth=40, seed=7)
    again = bt.pp.rarefy(bt.datasets.toy(), depth=40, seed=np.random.default_rng(7))
    assert (first.X != again.X).nnz == 0


def test_features_left_all_zero_are_dropped():
    dense = np.array([[5, 1, 0], [6, 0, 1]])
    out = bt.pp.rarefy(_adata(dense), depth=5, seed=0)
    assert out.X.toarray().sum(axis=0).min() > 0


def test_toy_tree_is_pruned_to_kept_features():
    out = bt.pp.rarefy(bt.datasets.toy(), depth=5, seed=0)
    tree = get_tree(out)
    assert {n for n in tree.nodes if tree.out_degree(n) == 0} == set(out.var_names)


def test_all_zero_sample_is_dropped_by_default():
    with pytest.warns(UserWarning, match=r"\['s0'\]"):
        out = bt.pp.rarefy(_adata(np.array([[0, 0], [3, 4], [5, 5]])), seed=0)
    assert list(out.obs_names) == ["s1", "s2"] and _depths(out) == [7, 7]


def test_single_sample():
    assert _depths(bt.pp.rarefy(_adata(np.array([[3, 9, 2]])), depth=10, seed=0)) == [10]


def test_drops_derived_slots_keeps_counts_and_records_provenance():
    out = bt.pp.rarefy(bt.pp.relative(bt.datasets.toy()), depth=60, seed=0)
    assert "relative" not in out.layers and out.uns["biotapy"]["x_kind"] == "counts"
    entry = json.loads(out.uns["biotapy"]["provenance"][-1])
    assert entry["step"] == "pp.rarefy" and entry["params"] == {"depth": 60}


def test_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.pp.rarefy(tdata, depth=60, seed=0)
    assert_unchanged(before, tdata)


def test_rejects_non_counts():
    tdata = bt.datasets.toy()
    tdata.uns["biotapy"]["x_kind"] = "relative"
    with pytest.raises(ValueError, match="pp.rarefy needs raw counts"):
        bt.pp.rarefy(tdata)


def test_rejects_fractional_values_labelled_counts():
    # subsample_counts would truncate them: [1.5, 2.5, 3.0] rarefied to 3 gives [0, 0, 3].
    tdata = bt.datasets.toy()
    tdata.X = sp.csr_matrix(tdata.X.toarray() + 0.5)
    with pytest.raises(ValueError, match="pp.rarefy needs raw counts in X, but X holds non-integer values"):
        bt.pp.rarefy(tdata, seed=0)


def test_depth_below_one_raises():
    with pytest.raises(ValueError, match="depth must be at least 1"):
        bt.pp.rarefy(bt.datasets.toy(), depth=0)


@pytest.mark.parametrize("depth", [5.5, True])
def test_non_integer_depth_raises(depth):
    with pytest.raises(TypeError, match="depth= must be an integer, got"):
        bt.pp.rarefy(bt.datasets.toy(), depth=depth)


def test_depth_above_every_sample_raises():
    with pytest.raises(ValueError, match="depth=1000"):
        bt.pp.rarefy(bt.datasets.toy(), depth=1000)


@given(
    arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 30)), st.integers(1, 40)
)
def test_kept_rows_sum_to_depth_and_never_exceed_the_original(dense, depth):
    deep = dense.sum(axis=1) >= depth
    if not deep.any():
        with pytest.raises(ValueError, match="depth="):
            bt.pp.rarefy(_adata(dense), depth=depth, seed=0)
        return
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # shallow samples are dropped with a warning
        out = bt.pp.rarefy(_adata(dense), depth=depth, seed=0)
    assert list(out.obs_names) == [f"s{i}" for i in np.flatnonzero(deep)]
    assert _depths(out) == [depth] * int(deep.sum())
    kept = dense[deep][:, [int(name[1:]) for name in out.var_names]]
    assert (out.X.toarray() <= kept).all()
