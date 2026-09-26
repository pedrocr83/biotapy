import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt
from biotapy._core import get_tree

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus"]
# Small taxonomy used only by test_matches_naive_phyloseq_reference (F1); kept separate
# from RANKS above, which mirrors the toy() dataset's six real rank columns.
_LINEAGE_RANKS = ("kingdom", "phylum", "class", "genus")
_TAXA_VALUES = st.sampled_from(["x", "y", None])


def _leaves(tdata) -> set[str]:
    tree = get_tree(tdata)
    return {n for n in tree.nodes if tree.out_degree(n) == 0}


def _naive_reference(x, var, rank, *, dropna):
    """phyloseq-style reference: group by lineage up to rank, keep the most abundant archetype."""
    upto = _LINEAGE_RANKS[: _LINEAGE_RANKS.index(rank) + 1]
    keep_mask = var[rank].notna() if dropna else pd.Series(True, index=var.index)
    names = list(var.index[keep_mask])
    kept = x[:, keep_mask.to_numpy()]
    lineage = ["|".join("NA" if pd.isna(v) else str(v) for v in var.loc[n, upto]) for n in names]
    totals = kept.sum(axis=0)
    columns = {}
    for key in dict.fromkeys(lineage):
        idx = [i for i, lin in enumerate(lineage) if lin == key]
        archetype = idx[int(np.argmax(totals[idx]))]
        columns[names[archetype]] = kept[:, idx].sum(axis=1)
    order = [n for n in names if n in columns]
    return order, np.column_stack([columns[n] for n in order])


@st.composite
def _tax_glom_cases(draw):
    n_obs = draw(st.integers(1, 4))
    n_var = draw(st.integers(1, 8))
    rank = draw(st.sampled_from(_LINEAGE_RANKS))
    counts = draw(st.lists(st.integers(0, 20), min_size=n_obs * n_var, max_size=n_obs * n_var))
    x = np.array(counts, dtype=np.int64).reshape(n_obs, n_var)
    var = pd.DataFrame(
        {r: draw(st.lists(_TAXA_VALUES, min_size=n_var, max_size=n_var)) for r in _LINEAGE_RANKS},
        index=[f"f{i}" for i in range(n_var)],
    )
    # dropna=True with no value anywhere at rank raises ValueError; that path has its own test.
    dropna = draw(st.booleans()) if var[rank].notna().any() else False
    return x, var, rank, dropna


def test_phylum_keeps_most_abundant_member_in_original_order():
    out = bt.pp.tax_glom(bt.datasets.toy(), "phylum")
    assert list(out.var_names) == ["f3", "f6", "f7"]
    assert out.var[["class", "order", "family", "genus"]].isna().all().all()
    assert _leaves(out) == {"f3", "f6", "f7"}


def test_genus_sums_members_and_drops_unassigned():
    out = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    assert list(out.var_names) == ["f1", "f2", "f3", "f4", "f6", "f7"]
    np.testing.assert_array_equal(out[:, "f4"].X.toarray().ravel(), [30, 25, 36, 10, 14, 9])


def test_dropna_false_keeps_unassigned_as_own_group():
    assert bt.pp.tax_glom(bt.datasets.toy(), "genus", dropna=False).n_vars == 7


def test_groups_by_lineage_not_label():
    tdata = bt.datasets.toy()
    tdata.var.loc["f7", "genus"] = "Blautia"
    assert bt.pp.tax_glom(tdata, "genus").n_vars == 6


def test_interleaved_taxa_columns_end_up_in_original_order():
    # a1/a2 share phylum A, b1 is phylum B; a2 is the most abundant member of A.
    # Column order must follow the archetypes' original positions (b1, a2), not group order.
    var = pd.DataFrame({"kingdom": ["K", "K", "K"], "phylum": ["A", "B", "A"]}, index=["a1", "b1", "a2"])
    x = np.array([[1, 5, 9], [2, 6, 8]])
    adata = AnnData(X=sp.csr_matrix(x), var=var, obs=pd.DataFrame(index=["s1", "s2"]))
    out = bt.pp.tax_glom(adata, "phylum")
    assert list(out.var_names) == ["b1", "a2"]
    np.testing.assert_array_equal(out.X.toarray(), [[5, 10], [6, 10]])


@given(_tax_glom_cases())
def test_matches_naive_phyloseq_reference(case):
    x, var, rank, dropna = case
    adata = AnnData(X=sp.csr_matrix(x), var=var, obs=pd.DataFrame(index=[f"s{i}" for i in range(x.shape[0])]))
    out = bt.pp.tax_glom(adata, rank, dropna=dropna)
    names, expected = _naive_reference(x, var, rank, dropna=dropna)
    assert list(out.var_names) == names
    np.testing.assert_array_equal(out.X.toarray(), expected)


def test_drops_derived_slots():
    assert "relative" not in bt.pp.tax_glom(bt.pp.relative(bt.datasets.toy()), "phylum").layers


def test_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.pp.tax_glom(tdata, "genus")
    assert_unchanged(before, tdata)


def test_unknown_rank_is_named():
    with pytest.raises(KeyError, match="species"):
        bt.pp.tax_glom(bt.datasets.toy(), "species")


def test_rank_with_no_values_raises():
    tdata = bt.datasets.toy()
    tdata.var["genus"] = None
    with pytest.raises(ValueError, match="genus"):
        bt.pp.tax_glom(tdata, "genus")


@given(st.sampled_from(RANKS))
def test_sample_totals_preserved_without_dropna(rank):
    tdata = bt.datasets.toy()
    out = bt.pp.tax_glom(tdata, rank, dropna=False)
    np.testing.assert_array_equal(np.asarray(out.X.sum(axis=1)).ravel(), np.asarray(tdata.X.sum(axis=1)).ravel())


def test_all_zero_sample_stays_zero():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[0] = 0
    tdata.X = sp.csr_matrix(dense)
    out = bt.pp.tax_glom(tdata, "phylum")
    assert out.X[0].nnz == 0
    unmodified = bt.pp.tax_glom(bt.datasets.toy(), "phylum").X
    np.testing.assert_array_equal(out.X[1:].toarray(), unmodified[1:].toarray())


def test_all_zero_group_keeps_first_member():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, [6, 7]] = 0  # f7, f8
    tdata.X = sp.csr_matrix(dense)
    out = bt.pp.tax_glom(tdata, "phylum")
    assert list(out.var_names) == ["f3", "f6", "f7"]
    np.testing.assert_array_equal(out[:, "f7"].X.toarray().ravel(), np.zeros(6))


def test_single_sample():
    tdata = bt.datasets.toy()[:1].copy()
    out = bt.pp.tax_glom(tdata, "phylum")
    assert list(out.var_names) == ["f3", "f4", "f7"]
    np.testing.assert_array_equal(out.X.toarray(), [[35, 32, 1]])
