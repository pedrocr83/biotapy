import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import skbio
from anndata import AnnData
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from mudata import MuData

import biotapy as bt


def _metabolites(counts: np.ndarray, samples) -> AnnData:
    """Three metabolites: one follows f6 (Prevotella), one f3 (Faecalibacterium), one is flat."""
    values = np.column_stack([counts[:, 5] + 1, counts[:, 2] + 1, np.full(counts.shape[0], 10)])
    return AnnData(
        X=sp.csr_matrix(values.astype(np.float64)),
        obs=pd.DataFrame(index=list(samples)),
        var=pd.DataFrame(index=["m_prev", "m_faec", "m_flat"]),
    )


def _toy_mudata():
    tdata = bt.datasets.toy()
    return bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata.X.toarray(), tdata.obs_names)})


def test_ranks_are_microbes_by_metabolites():
    ranks = bt.tl.mmvec(_toy_mudata(), seed=0)
    assert isinstance(ranks, pd.DataFrame)
    assert ranks.index.tolist() == [f"f{i}" for i in range(1, 9)]
    assert ranks.columns.tolist() == ["m_prev", "m_faec", "m_flat"]
    np.testing.assert_allclose(ranks.sum(axis=1), 0.0, atol=1e-10)


def test_a_metabolite_ranks_highest_for_the_microbe_it_follows():
    ranks = bt.tl.mmvec(_toy_mudata(), seed=0)
    assert ranks.loc["f6"].idxmax() == "m_prev"
    assert ranks.loc["f3"].idxmax() == "m_faec"


def test_same_seed_same_ranks():
    first = bt.tl.mmvec(_toy_mudata(), seed=1)
    again = bt.tl.mmvec(_toy_mudata(), seed=np.random.default_rng(1))
    pd.testing.assert_frame_equal(first, again)


def test_modality_names_are_arguments():
    mdata = _toy_mudata()
    renamed = bt.io.to_mudata({"microbes": mdata["taxa"], "compounds": mdata["metabolites"]})
    pd.testing.assert_frame_equal(
        bt.tl.mmvec(renamed, microbes="microbes", metabolites="compounds", seed=0), bt.tl.mmvec(mdata, seed=0)
    )


def test_ranks_stay_a_labelled_dataframe_whatever_scikit_bio_outputs():
    previous = skbio.get_config("table_output")
    skbio.set_config("table_output", "numpy")
    try:
        ranks = bt.tl.mmvec(_toy_mudata(), seed=0)
    finally:
        skbio.set_config("table_output", previous)
    assert isinstance(ranks, pd.DataFrame)
    assert ranks.columns.tolist() == ["m_prev", "m_faec", "m_flat"]


def test_keeps_the_input(assert_unchanged):
    mdata = _toy_mudata()
    before = {name: mod.copy() for name, mod in mdata.mod.items()}
    bt.tl.mmvec(mdata, seed=0)
    for name, mod in mdata.mod.items():
        assert_unchanged(before[name], mod)


def test_missing_modality_raises():
    with pytest.raises(KeyError, match=r"metabolites='metabolome' is not a modality; found \['taxa', 'metabolites'\]"):
        bt.tl.mmvec(_toy_mudata(), metabolites="metabolome")


def test_unaligned_samples_raise():
    tdata = bt.datasets.toy()
    metabolites = _metabolites(tdata.X.toarray(), tdata.obs_names)[::-1].copy()
    with pytest.raises(ValueError, match=r"hold different samples.*bt.io.to_mudata"):
        bt.tl.mmvec(MuData({"taxa": tdata, "metabolites": metabolites}))


def test_all_zero_feature_raises():
    mdata = _toy_mudata()
    mdata["taxa"].X = sp.csr_matrix(np.column_stack([mdata["taxa"].X.toarray()[:, :7], np.zeros(6)]))
    with pytest.raises(ValueError, match=r"microbes='taxa' has all-zero features: \['f8'\]"):
        bt.tl.mmvec(mdata)


def test_all_zero_sample_raises():
    mdata = _toy_mudata()
    values = mdata["metabolites"].X.toarray()
    values[2] = 0.0
    mdata["metabolites"].X = sp.csr_matrix(values)
    with pytest.raises(ValueError, match=r"metabolites='metabolites' has all-zero samples: \['s3'\]"):
        bt.tl.mmvec(mdata)


def test_single_sample():
    tdata = bt.pp.filter_features(bt.datasets.toy()[:1].copy(), min_total=1)
    metabolites = AnnData(
        X=sp.csr_matrix([[3.0, 4.0, 5.0]]), obs=pd.DataFrame(index=["s1"]), var=pd.DataFrame(index=["m1", "m2", "m3"])
    )
    assert bt.tl.mmvec(bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites}), seed=0).shape == (6, 3)


@pytest.mark.parametrize("bad", [-1.0, np.nan])
def test_negative_or_missing_value_raises(bad):
    mdata = _toy_mudata()
    values = mdata["metabolites"].X.toarray()
    values[0, 0] = bad
    mdata["metabolites"].X = sp.csr_matrix(values)
    with pytest.raises(ValueError, match="tl.mmvec needs finite, non-negative values in metabolites='metabolites'"):
        bt.tl.mmvec(mdata)


@settings(max_examples=25, deadline=None)
@given(arrays(np.int64, st.tuples(st.integers(2, 5), st.integers(2, 4)), elements=st.integers(1, 50)))
def test_every_microbe_row_is_centred(counts):
    samples = [f"s{i}" for i in range(counts.shape[0])]
    taxa = AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=samples),
        var=pd.DataFrame(index=[f"f{i}" for i in range(counts.shape[1])]),
    )
    metabolites = AnnData(
        X=sp.csr_matrix(counts[:, ::-1] + 1),
        obs=pd.DataFrame(index=samples),
        var=pd.DataFrame(index=[f"m{i}" for i in range(counts.shape[1])]),
    )
    ranks = bt.tl.mmvec(bt.io.to_mudata({"taxa": taxa, "metabolites": metabolites}), seed=0)
    np.testing.assert_allclose(ranks.sum(axis=1), 0.0, atol=1e-9)
