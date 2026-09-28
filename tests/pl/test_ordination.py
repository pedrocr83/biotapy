import numpy as np
import pytest

import biotapy as bt


def _toy_with(*bases):
    tdata = bt.datasets.toy()
    bt.tl.beta(tdata, inplace=True)
    if "pcoa" in bases:
        bt.tl.pcoa(tdata, inplace=True)
    if "nmds" in bases:
        bt.tl.nmds(tdata, seed=0, inplace=True)
    return tdata


def _points(ax):
    return np.concatenate([collection.get_offsets() for collection in ax.collections])


def test_ordination_draws_stored_pcoa_axes_with_percentages(ax):
    tdata = _toy_with("pcoa")
    bt.pl.ordination(tdata, ax=ax)
    np.testing.assert_array_equal(_points(ax), tdata.obsm["X_pcoa"][:, :2])
    share = tdata.uns["biotapy"]["pcoa"]["proportion_explained"]
    assert ax.get_xlabel() == f"PC1 [{100 * share[0]:.1f}%]"
    assert ax.get_ylabel() == f"PC2 [{100 * share[1]:.1f}%]"


def test_ordination_components_pick_the_axes(ax):
    tdata = _toy_with("pcoa")
    bt.pl.ordination(tdata, components=(1, 3), ax=ax)
    np.testing.assert_array_equal(_points(ax), tdata.obsm["X_pcoa"][:, [0, 2]])
    assert ax.get_ylabel().startswith("PC3 [")


def test_ordination_nmds_has_no_percentages_and_notes_the_stress(ax):
    tdata = _toy_with("nmds")
    bt.pl.ordination(tdata, basis="nmds", color="group", ax=ax)
    np.testing.assert_array_equal(np.sort(_points(ax), axis=0), np.sort(tdata.obsm["X_nmds"], axis=0))
    assert (ax.get_xlabel(), ax.get_ylabel()) == ("NMDS1", "NMDS2")
    assert ax.texts[0].get_text() == f"stress {tdata.uns['biotapy']['nmds']['stress']:.3f}"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["A", "B"]


def test_ordination_input_unchanged(assert_unchanged, ax):
    tdata = _toy_with("pcoa", "nmds")
    before = tdata.copy()
    bt.pl.ordination(tdata, basis="nmds", color="group", ax=ax)
    bt.pl.scree(tdata, ax=ax)
    assert_unchanged(before, tdata)


def test_ordination_missing_names_the_call(ax):
    with pytest.raises(KeyError, match=r"run bt.tl.pcoa\(adata, inplace=True\) first"):
        bt.pl.ordination(_toy_with(), ax=ax)
    with pytest.raises(KeyError, match=r"run bt.tl.nmds\(adata, inplace=True\) first"):
        bt.pl.ordination(_toy_with("pcoa"), basis="nmds", ax=ax)


def test_ordination_is_gone_after_a_feature_change(ax):
    tdata = bt.pp.filter_features(_toy_with("pcoa"), min_total=19)
    with pytest.raises(KeyError, match="bt.tl.pcoa"):
        bt.pl.ordination(tdata, ax=ax)


def test_ordination_bad_basis_or_components_raise(ax):
    tdata = _toy_with("pcoa")
    with pytest.raises(ValueError, match="basis must be one of"):
        bt.pl.ordination(tdata, basis="tsne", ax=ax)
    with pytest.raises(ValueError, match="components must be two axis numbers from 1 to 5"):
        bt.pl.ordination(tdata, components=(1, 6), ax=ax)
    with pytest.raises(ValueError, match="components"):
        bt.pl.ordination(tdata, components=(0, 1), ax=ax)


def test_scree_draws_the_stored_proportions(ax):
    tdata = _toy_with("pcoa")
    bt.pl.scree(tdata, ax=ax)
    heights = [patch.get_height() for patch in ax.patches]
    np.testing.assert_array_equal(heights, tdata.uns["biotapy"]["pcoa"]["proportion_explained"])
    assert [label.get_text() for label in ax.get_xticklabels()] == ["PC1", "PC2", "PC3", "PC4", "PC5"]


def test_scree_missing_pcoa_names_the_call(ax):
    with pytest.raises(KeyError, match="bt.tl.pcoa"):
        bt.pl.scree(_toy_with("nmds"), ax=ax)
