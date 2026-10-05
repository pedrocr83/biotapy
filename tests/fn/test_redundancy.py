import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt

COLUMNS = ["taxonomic_diversity", "functional_diversity", "redundancy", "normalized_redundancy"]
# Three genomes over three genes. Weighted Jaccard (Tian et al. Eq. 7): d(A,B) = 1 - 2/4 = 0.5,
# d(A,C) = 1 - 0/4 = 1, d(B,C) = 1 - 1/5 = 0.8.
TRAITS = pd.DataFrame([[2, 1, 0], [1, 1, 1], [0, 0, 3]], index=["A", "B", "C"], columns=["g1", "g2", "g3"])


def _adata(dense, taxa=("A", "B", "C")):
    dense = np.asarray(dense, dtype=np.float64)
    return ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=list(taxa)),
    )


def test_hand_computed_samples():
    out = bt.fn.functional_redundancy(_adata([[2, 1, 1], [0, 3, 1]]), traits=TRAITS)
    assert out.columns.tolist() == COLUMNS and out.index.tolist() == ["s0", "s1"]
    # s0: p = (1/2, 1/4, 1/4). TD = 1 - 3/8 = 0.625; FD = 2 (0.5/8 + 1/8 + 0.8/16) = 0.475; FR = 0.15.
    np.testing.assert_allclose(out.loc["s0"].to_numpy(), [0.625, 0.475, 0.15, 0.24], rtol=1e-12)
    # s1: p = (0, 3/4, 1/4). TD = 0.375; FD = 2 x 0.8 x 3/16 = 0.3; FR = 0.075; nFR = 0.2.
    np.testing.assert_allclose(out.loc["s1"].to_numpy(), [0.375, 0.3, 0.075, 0.2], rtol=1e-12)


def test_taxonomic_diversity_is_gini_simpson():
    adata = _adata([[2, 1, 1], [0, 3, 1], [5, 0, 1]])
    out = bt.fn.functional_redundancy(adata, traits=TRAITS)
    simpson = bt.tl.alpha(adata, metrics=["simpson"])["simpson"]
    np.testing.assert_allclose(out["taxonomic_diversity"].to_numpy(), simpson.to_numpy(), rtol=1e-12)


def test_scale_of_each_sample_does_not_matter():
    out = bt.fn.functional_redundancy(_adata([[2, 1, 1], [200, 100, 100]]), traits=TRAITS)
    np.testing.assert_allclose(out.loc["s0"].to_numpy(), out.loc["s1"].to_numpy(), rtol=1e-12)


def test_taxa_without_traits_are_left_out_with_a_warning():
    adata = _adata([[2, 1, 1, 7], [0, 3, 1, 0]], taxa=("A", "B", "C", "RARE"))
    with pytest.warns(UserWarning, match=r"1 of 4 taxa have no row in traits.*\['RARE'\]"):
        out = bt.fn.functional_redundancy(adata, traits=TRAITS)
    expected = bt.fn.functional_redundancy(_adata([[2, 1, 1], [0, 3, 1]]), traits=TRAITS)
    pd.testing.assert_frame_equal(out, expected)


def test_extra_trait_rows_and_order_do_not_matter():
    traits = pd.concat([TRAITS, pd.DataFrame([[9, 9, 9]], index=["Z"], columns=TRAITS.columns)]).iloc[::-1]
    out = bt.fn.functional_redundancy(_adata([[2, 1, 1]]), traits=traits)
    np.testing.assert_allclose(out.loc["s0"].to_numpy(), [0.625, 0.475, 0.15, 0.24], rtol=1e-12)


def test_integer_trait_ids_match_text_var_names():
    traits = TRAITS.set_axis([1, 2, 3])
    out = bt.fn.functional_redundancy(_adata([[2, 1, 1]], taxa=("1", "2", "3")), traits=traits)
    np.testing.assert_allclose(out.loc["s0", "redundancy"], 0.15, rtol=1e-12)


def test_reads_picrust2_traits_and_16s_corrected_abundances(tmp_path):
    (tmp_path / "EC_predicted.tsv").write_text(
        "sequence\tEC:1.1.1.1\tEC:2.7.1.1\tEC:3.2.1.1\tmetadata_NSTI\nASV1\t2\t1\t0\t0.1\n0042\t1\t1\t1\t0.2\n"
    )
    (tmp_path / "marker.tsv").write_text("sequence\t16S_rRNA_Count\tmetadata_NSTI\nASV1\t2\t0.1\n0042\t1\t0.2\n")
    traits = bt.io.read_picrust2_traits(tmp_path / "EC_predicted.tsv")
    copies = bt.io.read_picrust2_traits(tmp_path / "marker.tsv")["16S_rRNA_Count"]
    reads = _adata([[4, 2]], taxa=("ASV1", "0042"))
    cells = reads.copy()
    cells.X = sp.csr_matrix(reads.X.multiply(1 / copies[reads.var_names].to_numpy()))
    # Cells: ASV1 4/2 = 2, 0042 2/1 = 2, so p = (1/2, 1/2); d = 0.5, TD = 0.5, FD = 0.25.
    out = bt.fn.functional_redundancy(cells, traits=traits)
    np.testing.assert_allclose(out.loc["s0"].to_numpy(), [0.5, 0.25, 0.25, 0.5], rtol=1e-12)


def test_inputs_unchanged(assert_unchanged):
    adata, traits = _adata([[2, 1, 1], [0, 3, 1]]), TRAITS.copy()
    before = adata.copy()
    bt.fn.functional_redundancy(adata, traits=traits)
    assert_unchanged(before, adata)
    pd.testing.assert_frame_equal(traits, TRAITS)


def test_all_zero_sample_is_nan():
    out = bt.fn.functional_redundancy(_adata([[0, 0, 0], [2, 1, 1]]), traits=TRAITS)
    assert out.loc["s0"].isna().all() and out.loc["s1"].notna().all()


def test_sample_whose_taxa_all_lack_traits_is_nan():
    adata = _adata([[0, 0, 0, 5], [2, 1, 1, 0]], taxa=("A", "B", "C", "RARE"))
    with pytest.warns(UserWarning, match="no row in traits"):
        out = bt.fn.functional_redundancy(adata, traits=TRAITS)
    assert out.loc["s0"].isna().all()


def test_single_taxon_sample_has_no_diversity_and_undefined_normalized_redundancy():
    out = bt.fn.functional_redundancy(_adata([[0, 0, 5]]), traits=TRAITS)
    assert out.loc["s0", COLUMNS[:3]].tolist() == [0.0, 0.0, 0.0] and np.isnan(out.loc["s0", COLUMNS[3]])


def test_all_zero_taxon_changes_nothing():
    out = bt.fn.functional_redundancy(
        _adata([[2, 1, 1, 0]], taxa=("A", "B", "C", "D")), traits=TRAITS.reindex(["A", "B", "C", "D"], fill_value=1)
    )
    np.testing.assert_allclose(out.loc["s0"].to_numpy(), [0.625, 0.475, 0.15, 0.24], rtol=1e-12)


def test_taxa_without_genes_share_nothing():
    traits = pd.DataFrame(np.zeros((2, 3)), index=["E", "F"], columns=["g1", "g2", "g3"])
    out = bt.fn.functional_redundancy(_adata([[1, 1]], taxa=("E", "F")), traits=traits)
    assert out.loc["s0", COLUMNS].tolist() == [0.5, 0.5, 0.0, 0.0]


def test_single_sample():
    assert bt.fn.functional_redundancy(_adata([[2, 1, 1]]), traits=TRAITS).shape == (1, 4)


def test_no_shared_ids_raise_with_examples():
    with pytest.raises(ValueError, match=r"adata: \['f1', 'f2', 'f3'\], traits: \['A', 'B', 'C'\]"):
        bt.fn.functional_redundancy(_adata([[2, 1, 1]], taxa=("f1", "f2", "f3")), traits=TRAITS)


@pytest.mark.parametrize(
    ("traits", "error", "message"),
    [
        (TRAITS.to_numpy(), TypeError, "traits must be a pandas.DataFrame"),
        (TRAITS.assign(g4="x"), TypeError, r"non-numeric columns: \['g4'\]"),
        (TRAITS.set_axis(["A", "A", "C"]), ValueError, r"traits repeats row ids: \['A'\]"),
        (TRAITS.assign(g1=[np.nan, 1, 0]), ValueError, "traits holds a missing, negative or infinite"),
        (TRAITS.assign(g1=[-1, 1, 0]), ValueError, "traits holds a missing, negative or infinite"),
        (TRAITS.assign(g1=[np.inf, 1, 0]), ValueError, "traits holds a missing, negative or infinite"),
    ],
    ids=["array", "text-column", "repeated-id", "nan", "negative", "infinite"],
)
def test_bad_traits_raise_naming_traits(traits, error, message):
    with pytest.raises(error, match=message):
        bt.fn.functional_redundancy(_adata([[2, 1, 1]]), traits=traits)


@pytest.mark.parametrize("value", [-1.0, np.nan, np.inf])
def test_bad_abundances_raise_naming_adata(value):
    with pytest.raises(ValueError, match="adata: X holds a missing, negative or infinite"):
        bt.fn.functional_redundancy(_adata([[2, 1, value]]), traits=TRAITS)


def test_overflowing_sample_total_raises_naming_adata():
    with pytest.raises(ValueError, match=r"adata: the total abundance of a sample overflows: \['s0'\]"):
        bt.fn.functional_redundancy(_adata([[1e308, 1e308, 0], [1, 1, 1]]), traits=TRAITS)


def test_repeated_taxon_ids_raise_naming_adata():
    with pytest.warns(UserWarning, match="not unique"):  # AnnData itself warns on construction
        adata = _adata([[2, 1, 1]], taxa=("A", "A", "C"))
    with pytest.raises(ValueError, match=r"adata repeats taxon ids: \['A'\]"):
        bt.fn.functional_redundancy(adata, traits=TRAITS)


@pytest.mark.parametrize("dtype", ["Int64", "Float64"])
def test_nullable_missing_trait_raises_naming_traits(dtype):
    traits = TRAITS.astype(dtype)
    traits.iloc[0, 0] = pd.NA
    with pytest.raises(ValueError, match="traits holds a missing, negative or infinite"):
        bt.fn.functional_redundancy(_adata([[2, 1, 1]]), traits=traits)


@st.composite
def _cases(draw):
    n_taxa, n_genes, n_obs = draw(st.integers(1, 6)), draw(st.integers(1, 5)), draw(st.integers(1, 3))
    genomes = draw(arrays(np.int64, (n_taxa, n_genes), elements=st.integers(0, 4)))
    abundances = draw(arrays(np.float64, (n_obs, n_taxa), elements=st.floats(0, 1e3)))
    return genomes, abundances


def _run(genomes, abundances, order=None):
    order = np.arange(genomes.shape[0]) if order is None else order
    taxa = [f"t{i}" for i in order]
    traits = pd.DataFrame(genomes[order], index=taxa)
    return bt.fn.functional_redundancy(_adata(abundances[:, order], taxa=taxa), traits=traits)


def _weighted_jaccard(u, v):
    # Eq. 7 written out, independent of SciPy; two empty genomes share nothing.
    larger = np.maximum(u, v).sum()
    return 1 - np.minimum(u, v).sum() / larger if larger > 0 else 1.0


@settings(deadline=None)
@given(_cases())
def test_diversities_are_ordered(case):
    out = _run(*case).dropna(subset=["taxonomic_diversity"])
    td, fd, fr = (out[column].to_numpy() for column in COLUMNS[:3])
    assert (fd >= 0).all() and (fr >= 0).all() and (fd <= td).all() and (fr <= td).all()
    np.testing.assert_allclose(fd + fr, td, rtol=1e-12, atol=1e-15)


@settings(deadline=None)
@given(_cases())
def test_matches_the_paper_equations(case):
    genomes, abundances = case
    out = _run(genomes, abundances)
    for sample, row in enumerate(abundances):
        if row.sum() == 0:
            assert out.iloc[sample].isna().all()
            continue
        p = row / row.sum()
        pairs = [(i, j) for i in range(p.size) for j in range(p.size) if i != j]
        fd = sum(_weighted_jaccard(genomes[i], genomes[j]) * p[i] * p[j] for i, j in pairs)
        # atol: values are sums of products of shares, so near-zero results carry absolute rounding.
        np.testing.assert_allclose(out.iloc[sample, :3], [1 - p @ p, fd, 1 - p @ p - fd], rtol=1e-9, atol=1e-12)


@settings(deadline=None)
@given(_cases(), st.randoms(use_true_random=False))
def test_taxon_order_does_not_matter(case, random):
    genomes, abundances = case
    order = np.array(random.sample(range(genomes.shape[0]), genomes.shape[0]))
    pd.testing.assert_frame_equal(_run(genomes, abundances, order), _run(genomes, abundances), rtol=1e-12)


@settings(deadline=None)
@given(_cases())
def test_identical_genomes_are_fully_redundant(case):
    genomes, abundances = case
    out = _run(np.tile(genomes[:1] + 1, (genomes.shape[0], 1)), abundances)
    np.testing.assert_array_equal(out["redundancy"].to_numpy(), out["taxonomic_diversity"].to_numpy())


@settings(deadline=None)
@given(_cases())
def test_disjoint_genomes_have_no_redundancy(case):
    genomes, abundances = case
    out = _run(np.diag(np.arange(1, genomes.shape[0] + 1)), abundances)
    assert (out["redundancy"].dropna() == 0).all()


def test_mudata_adata_raises_naming_adata_and_the_modality():
    import mudata

    mdata = mudata.MuData({"function": _adata([[2, 1, 1]])})
    with pytest.raises(TypeError, match=r"adata must be an AnnData.*MuData.*pass one modality"):
        bt.fn.functional_redundancy(mdata, traits=TRAITS)


def test_dataframe_adata_raises_naming_adata():
    with pytest.raises(TypeError, match=r"adata must be an AnnData.*DataFrame"):
        bt.fn.functional_redundancy(pd.DataFrame([[2, 1, 1]], columns=["A", "B", "C"]), traits=TRAITS)  # type: ignore[arg-type]
