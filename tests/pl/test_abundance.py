import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from matplotlib.colors import to_hex
from matplotlib.figure import Figure

import biotapy as bt
from biotapy._core import make_function_mudata


def _adata(dense, var=None) -> ad.AnnData:
    return ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(var, index=[f"f{i}" for i in range(dense.shape[1])]),
    )


def _heights(ax, n_bars):
    # ax.bar adds one Rectangle per bar, one group after another: reshape to bars x groups.
    return np.array([patch.get_height() for patch in ax.patches]).reshape(-1, n_bars).T


def _legend(ax):
    return [text.get_text() for text in ax.get_legend().get_texts()]


def _by_phylum(tdata):
    frame = pd.DataFrame(tdata.X.toarray(), columns=tdata.var["phylum"].to_numpy())
    return frame.T.groupby(level=0).sum().T  # samples x phyla, phyla sorted


def test_bar_sums_each_fill_group_per_sample(ax):
    tdata = bt.datasets.toy()
    bt.pl.bar(tdata, "phylum", ax=ax)
    expected = _by_phylum(tdata)
    assert _legend(ax) == list(expected.columns) == ["Bacteroidota", "Firmicutes", "Proteobacteria"]
    np.testing.assert_array_equal(_heights(ax, 6), expected.to_numpy())
    bottoms = np.array([patch.get_y() for patch in ax.patches]).reshape(-1, 6).T
    np.testing.assert_array_equal(bottoms[:, 1:], np.cumsum(expected.to_numpy(), axis=1)[:, :-1])
    assert [label.get_text() for label in ax.get_xticklabels()] == list(tdata.obs_names)


def test_bar_x_sums_the_samples_of_each_group(ax):
    tdata = bt.datasets.toy()
    bt.pl.bar(tdata, "phylum", x="group", ax=ax)
    expected = _by_phylum(tdata).groupby(tdata.obs["group"].to_numpy()).sum()
    np.testing.assert_array_equal(_heights(ax, 2), expected.to_numpy())
    assert [label.get_text() for label in ax.get_xticklabels()] == ["A", "B"]


def test_bar_fill_from_obs_colours_whole_samples(ax):
    tdata = bt.datasets.toy()
    bt.pl.bar(tdata, "group", ax=ax)
    totals = tdata.X.sum(axis=1).A1
    np.testing.assert_array_equal(_heights(ax, 6), np.c_[np.r_[totals[:3], 0, 0, 0], np.r_[0, 0, 0, totals[3:]]])
    assert _legend(ax) == ["A", "B"]


def test_bar_missing_rank_is_its_own_group_and_bars_keep_every_read(ax):
    tdata = bt.datasets.toy()  # f8 has no genus
    bt.pl.bar(tdata, "genus", ax=ax)
    assert _legend(ax)[-1] == "NA"
    np.testing.assert_array_equal(_heights(ax, 6)[:, -1], tdata.X[:, 7].toarray().ravel())
    np.testing.assert_array_equal(_heights(ax, 6).sum(axis=1), tdata.X.sum(axis=1).A1)


def test_bar_missing_obs_value_is_its_own_group(ax):
    tdata = bt.datasets.toy()
    tdata.obs["site"] = pd.Categorical(["u", None, "v", "u", "v", None])
    bt.pl.bar(tdata, "phylum", x="site", ax=ax)
    assert [label.get_text() for label in ax.get_xticklabels()] == ["u", "v", "NA"]
    np.testing.assert_array_equal(_heights(ax, 3).sum(), tdata.X.sum())


def test_bar_layer_plots_relative_abundance(ax):
    bt.pl.bar(bt.pp.relative(bt.datasets.toy()), "phylum", layer="relative", ax=ax)
    np.testing.assert_allclose(_heights(ax, 6).sum(axis=1), 1.0)
    assert ax.get_ylabel() == "relative"


def test_bar_all_zero_sample_and_feature(ax):
    dense = np.array([[0, 0, 0], [1, 0, 2], [3, 0, 0]])
    bt.pl.bar(_adata(dense, {"phylum": ["a", "b", "a"]}), "phylum", ax=ax)
    np.testing.assert_array_equal(_heights(ax, 3), [[0, 0], [3, 0], [3, 0]])


def test_bar_single_sample(ax):
    bt.pl.bar(_adata(np.array([[4, 0, 1]]), {"phylum": ["a", "b", "a"]}), "phylum", ax=ax)
    np.testing.assert_array_equal(_heights(ax, 1), [[5, 0]])


def test_bar_input_unchanged(assert_unchanged, ax):
    tdata = bt.pp.relative(bt.datasets.toy())
    before = tdata.copy()
    bt.pl.bar(tdata, "genus", x="group", layer="relative", ax=ax)
    assert_unchanged(before, tdata)


def test_bar_without_ax_draws_on_a_new_pyplot_figure():
    ax = bt.pl.bar(bt.datasets.toy(), "phylum")
    assert ax.figure.number in plt.get_fignums()
    plt.close(ax.figure)


def test_bar_missing_layer_names_the_call(ax):
    with pytest.raises(KeyError, match=r"run adata = bt\.pp\.relative\(adata\) first"):
        bt.pl.bar(bt.datasets.toy(), "phylum", layer="relative", ax=ax)


def test_bar_fill_must_be_in_exactly_one_of_var_and_obs(ax):
    with pytest.raises(KeyError, match="fill='site'"):
        bt.pl.bar(bt.datasets.toy(), "site", ax=ax)
    tdata = bt.datasets.toy()
    tdata.obs["phylum"] = "x"
    with pytest.raises(KeyError, match="exactly one of var and obs"):
        bt.pl.bar(tdata, "phylum", ax=ax)


def test_bar_numeric_x_raises(ax):
    tdata = bt.datasets.toy()
    tdata.obs["depth"] = tdata.X.sum(axis=1).A1
    with pytest.raises(TypeError, match=r"x='depth' is a numeric column"):
        bt.pl.bar(tdata, "phylum", x="depth", ax=ax)
    with pytest.raises(KeyError, match="x='site' is not a column of obs"):
        bt.pl.bar(tdata, "phylum", x="site", ax=ax)


@given(
    arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 20)),
    st.lists(st.sampled_from(["a", "b", None]), min_size=6, max_size=6),
)
def test_bar_heights_add_up_to_sample_totals(dense, labels):
    ax = Figure().add_subplot()
    bt.pl.bar(_adata(dense, {"phylum": labels[: dense.shape[1]]}), "phylum", ax=ax)
    np.testing.assert_array_equal(_heights(ax, dense.shape[0]).sum(axis=1), dense.sum(axis=1))


def test_heatmap_draws_the_table_features_by_samples(ax):
    tdata = bt.datasets.toy()
    image = bt.pl.heatmap(tdata, ax=ax).images[0]
    np.testing.assert_array_equal(image.get_array(), tdata.X.toarray().T)
    assert [label.get_text() for label in ax.get_xticklabels()] == list(tdata.obs_names)
    assert [label.get_text() for label in ax.get_yticklabels()] == list(tdata.var_names)
    assert len(ax.figure.axes) == 2  # the heatmap and its colour bar


def test_heatmap_zeros_are_drawn_black_on_a_log_scale(ax):
    image = bt.pl.heatmap(bt.datasets.toy(), ax=ax).images[0]
    assert np.ma.is_masked(image.norm(np.array([0.0])))  # LogNorm masks zeros; the colormap draws them "bad"
    assert image.cmap.get_bad().tolist() == [0.0, 0.0, 0.0, 1.0]


def test_heatmap_labels_at_most_250_names(ax):
    bt.pl.heatmap(_adata(np.ones((1, 251))), ax=ax)
    assert ax.get_yticks().size == 0
    at_cap = bt.pl.heatmap(_adata(np.ones((1, 250))), ax=Figure().add_subplot())
    assert [label.get_text() for label in at_cap.get_yticklabels()] == [f"f{i}" for i in range(250)]


def test_bar_labels_at_most_250_names(ax):
    bt.pl.bar(_adata(np.ones((251, 1)), var={"phylum": ["p"]}), "phylum", ax=ax)
    assert ax.get_xticks().size == 0
    at_cap = bt.pl.bar(_adata(np.ones((250, 1)), var={"phylum": ["p"]}), "phylum", ax=Figure().add_subplot())
    assert [label.get_text() for label in at_cap.get_xticklabels()] == [f"s{i}" for i in range(250)]


def test_heatmap_layer_and_all_zero_sample(ax):
    tdata = bt.pp.relative(bt.datasets.toy())
    image = bt.pl.heatmap(tdata, layer="relative", ax=ax).images[0]
    np.testing.assert_allclose(image.get_array().sum(axis=0), 1.0)
    dense = np.array([[0, 0], [1, 2]])
    np.testing.assert_array_equal(
        bt.pl.heatmap(_adata(dense), ax=Figure().add_subplot()).images[0].get_array(), dense.T
    )


def test_heatmap_single_sample(ax):
    dense = np.array([[4, 0, 1]])
    np.testing.assert_array_equal(bt.pl.heatmap(_adata(dense), ax=ax).images[0].get_array(), dense.T)


def test_heatmap_nothing_positive_raises(ax):
    with pytest.raises(ValueError, match=r"^adata: X holds no positive value"):
        bt.pl.heatmap(_adata(np.zeros((2, 3))), ax=ax)
    adata = _adata(np.ones((2, 3)))
    adata.layers["relative"] = sp.csr_matrix((2, 3))
    with pytest.raises(ValueError, match=r"^layer='relative': layers\['relative'\] holds no positive value"):
        bt.pl.heatmap(adata, layer="relative", ax=ax)


def test_heatmap_input_unchanged(assert_unchanged, ax):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.pl.heatmap(tdata, ax=ax)
    assert_unchanged(before, tdata)


def _by_taxon():
    return bt.datasets.toy_humann()["function_by_taxon"]


def test_contributions_stacks_one_segment_per_taxon(ax):
    by_taxon = _by_taxon()
    assert bt.pl.contributions(by_taxon, "2.7.1.2", ax=ax) is ax
    expected = bt.fn.contributions(by_taxon, "2.7.1.2")
    np.testing.assert_array_equal(_heights(ax, 6), expected.to_numpy())
    assert _legend(ax) == expected.columns.tolist() and ax.get_legend().get_title().get_text() == "taxon"
    assert [label.get_text() for label in ax.get_xticklabels()] == list(by_taxon.obs_names)
    assert ax.get_title() == "2.7.1.2" and ax.get_xlabel() == "sample"


def test_contributions_top_draws_other_last_in_grey(ax):
    bt.pl.contributions(_by_taxon(), "2.7.1.2", top=1, ax=ax)
    assert _legend(ax) == ["g__Blautia.s__Blautia_obeum", "other"]
    assert to_hex(ax.patches[-1].get_facecolor()) == "#7f7f7f"
    np.testing.assert_array_equal(_heights(ax, 6)[:, 1], [6, 5, 7, 1, 0, 2])


def test_contributions_a_taxon_named_other_keeps_its_colour(ax):
    ids = pd.Index(["K1", "K1|other", "K1|a"])
    mdata = make_function_mudata(
        np.array([[2, 1, 1]]), obs=pd.DataFrame(index=["s1"]), row_ids=ids, x_kind="rpk", source="t"
    )
    bt.pl.contributions(mdata["function_by_taxon"], "K1", ax=ax)
    assert "#7f7f7f" not in [to_hex(patch.get_facecolor()) for patch in ax.patches]


def test_contributions_draws_regrouped_tables(ax):
    edges = pd.DataFrame({"child": ["2.7.1.1", "2.7.1.2"], "parent": "kinase", "level": "role"})
    by_role = bt.fn.func_glom(_by_taxon(), "role", hierarchy=edges)
    bt.pl.contributions(by_role, "kinase", ax=ax)
    np.testing.assert_array_equal(_heights(ax, 6), bt.fn.contributions(by_role, "kinase").to_numpy())


def test_contributions_all_zero_sample_and_single_sample(ax):
    by_taxon = _by_taxon()
    dense = by_taxon.X.toarray()
    dense[0] = 0
    by_taxon.X = sp.csr_matrix(dense)
    bt.pl.contributions(by_taxon, "2.7.1.2", ax=ax)
    assert _heights(ax, 6)[0].tolist() == [0.0, 0.0]
    single = bt.pl.contributions(_by_taxon()[:1].copy(), "2.7.1.2", ax=Figure().add_subplot())
    assert len(single.patches) == 2


def test_contributions_input_unchanged(assert_unchanged, ax):
    by_taxon = _by_taxon()
    before = by_taxon.copy()
    bt.pl.contributions(by_taxon, "2.7.1.2", top=1, ax=ax)
    assert_unchanged(before, by_taxon)


def test_contributions_without_ax_draws_on_a_new_pyplot_figure():
    ax = bt.pl.contributions(_by_taxon(), "2.7.1.2")
    assert ax.figure.number in plt.get_fignums()
    plt.close(ax.figure)


def test_contributions_errors_name_the_argument(ax):
    with pytest.raises(KeyError, match="function='2.7.1.3'"):
        bt.pl.contributions(_by_taxon(), "2.7.1.3", ax=ax)
    with pytest.raises(ValueError, match="top=0"):
        bt.pl.contributions(_by_taxon(), "2.7.1.2", top=0, ax=ax)


def _function_with_taxa(n_taxa):
    ids = pd.Index(["K1", *[f"K1|t{j:02d}" for j in range(n_taxa)]])
    dense = np.arange(1, len(ids) + 1, dtype=np.float64)[None, :] * np.array([[1.0], [2.0]])
    mdata = make_function_mudata(dense, obs=pd.DataFrame(index=["s1", "s2"]), row_ids=ids, x_kind="rpk", source="t")
    return mdata["function_by_taxon"]


def test_contributions_other_is_not_the_colour_of_any_taxon(ax):
    bt.pl.contributions(_function_with_taxa(12), "K1", ax=ax)
    colours = [to_hex(patch.get_facecolor()) for patch in ax.patches[::2]]  # one patch per segment, sample s1 first
    assert len(colours) == 9
    assert colours[-1] == "#7f7f7f"
    assert "#7f7f7f" not in colours[:-1]
    assert len(set(colours)) == 9


def test_bar_na_group_has_a_colour_no_group_shares(ax):
    adata = _adata(np.ones((2, 9)), var={"g": [f"g{i}" for i in range(8)] + [None]})
    adata.var["g"] = pd.Categorical(adata.var["g"])
    bt.pl.bar(adata, "g", ax=ax)
    colours = [to_hex(patch.get_facecolor()) for patch in ax.patches[::2]]
    assert colours[-1] == "#7f7f7f"
    assert "#7f7f7f" not in colours[:-1]
    assert len(set(colours)) == 9


@pytest.mark.parametrize("n", [1, 5, 7])
def test_colours_with_na_keep_the_first_seven_of_tab10(n):
    from matplotlib import colormaps

    from biotapy.pl._common import _colors

    assert _colors(n, missing=True)[:n] == [colormaps["tab10"](i) for i in range(n)]


@pytest.mark.parametrize("n", [3, 10, 11, 20, 25])
def test_colours_without_na_are_tab10_tab20_or_turbo(n):
    from matplotlib import colormaps

    from biotapy.pl._common import _colors

    cmap = colormaps["tab10"] if n <= 10 else colormaps["tab20"] if n <= 20 else colormaps["turbo"].resampled(n)
    assert _colors(n, missing=False) == [cmap(i) for i in range(n)]


@pytest.mark.parametrize("n", [8, 9, 10, 15, 19, 20, 25])
def test_colours_with_na_never_hold_the_missing_grey(n):
    from biotapy.pl._common import _colors

    colours = [to_hex(colour) for colour in _colors(n, missing=True)]
    assert colours[-1] == "#7f7f7f"
    assert "#7f7f7f" not in colours[:-1]
    assert len(set(colours)) == n + 1


def test_contributions_a_mudata_or_the_community_modality_raises_as_fn_does(ax):
    mdata = bt.datasets.toy_humann()
    with pytest.raises(TypeError, match=r"mdata\['function_by_taxon'\]"):
        bt.pl.contributions(mdata, "2.7.1.2", ax=ax)
    with pytest.raises(KeyError, match=r"adata needs var columns \['function', 'taxon'\].*function_by_taxon"):
        bt.pl.contributions(mdata["function"], "2.7.1.2", ax=ax)
