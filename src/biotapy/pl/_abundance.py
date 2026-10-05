"""Abundance plots: stacked bars, a heatmap of the table, and the taxa behind one function."""

from typing import TYPE_CHECKING, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import sum_by
from biotapy.fn import contributions as function_contributions

from ._common import RGBA, _colors, groups, label_ticks, new_axes, obs_groups, table

if TYPE_CHECKING:
    from matplotlib.axes import Axes

# phyloseq's plot_heatmap colours: low "#000033", high "#66CCFF", and na.value "black", which zeros become on its log scale.
_LOW, _HIGH, _ZERO = "#000033", "#66CCFF", "black"


def bar(
    adata: AnnData, fill: str, *, x: str | None = None, layer: str | None = None, ax: "Axes | None" = None
) -> "Axes":
    """Stacked bars of abundance, one segment per group of ``fill``.

    Parameters
    ----------
    adata
        Samples x features.
    fill
        A ``var`` column, such as a taxonomic rank, or an ``obs`` column; each of its
        values is one coloured segment.
    x
        An ``obs`` column: one bar per value, summing its samples. By default one bar per sample.
    layer
        Plot ``layers[layer]``, such as ``"relative"`` from :func:`biotapy.pp.relative`, instead of ``X``.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        Bars in ``obs`` order (samples) or group order (``x``), segments in group order,
        and a legend titled ``fill``.

    Raises
    ------
    KeyError
        ``fill`` is in neither or both of ``var`` and ``obs``, ``x`` is not an ``obs``
        column, or ``layers[layer]`` is missing.
    TypeError
        ``fill`` or ``x`` is a numeric column.
    ValueError
        The plotted table holds a negative value, as ``layers["clr"]`` does.

    Notes
    -----
    R equivalent: ``phyloseq::plot_bar``
    Guide: :doc:`/guide/plotting`

    phyloseq stacks one outlined segment per feature and sample, which is 500,000
    rectangles on GlobalPatterns. biotapy sums the features of each ``fill`` group
    first: bar and segment heights are the same, without the per-feature outlines.
    Groups are categories in their order, else sorted values; missing values form
    the last group, ``NA``, so bars still sum every feature. ``facet_grid`` has no
    equivalent: draw one subset per axes of ``matplotlib.pyplot.subplots``.

    Examples
    --------
    >>> import biotapy as bt
    >>> ax = bt.pl.bar(bt.datasets.toy(), "phylum")
    >>> len(ax.patches)  # 6 samples x 3 phyla
    18
    """
    values = table(adata, layer, func="pl.bar")
    segments, labels, colors = _segments(adata, values, fill)
    heights, names = _by_x(adata, segments, x)
    ax = new_axes(ax)
    _stack(ax, heights, labels, colors=colors)
    label_ticks(ax, names, axis="x")
    ax.set_xlabel(x or "sample")
    ax.set_ylabel(layer or "abundance")
    ax.legend(title=fill)
    return ax


def contributions(adata: AnnData, function: str, *, top: int | None = 8, ax: "Axes | None" = None) -> "Axes":
    """Stacked bars of one function's abundance per sample, one segment per taxon.

    Parameters
    ----------
    adata
        A stratified function table: the ``"function_by_taxon"`` modality of
        ``bt.io.read_humann`` or ``bt.io.read_picrust2``, or
        ``bt.fn.func_glom``'s or ``bt.fn.renorm``'s output for it.
    function
        A function id as it appears in ``var["function"]``.
    top
        Draw the ``top`` taxa with the largest total over all samples and sum
        the rest into a last, grey segment, ``"other"``; ``None`` draws every
        taxon.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        One bar per sample in ``obs`` order, segments as the columns of
        ``bt.fn.contributions(adata, function, top=top)``, a legend titled
        ``taxon`` and the function id as the title.

    Raises
    ------
    TypeError
        ``adata`` is not an AnnData.
    KeyError
        ``adata`` lacks the ``function`` or ``taxon`` column; ``function`` has
        no stratified row (the message lists close ids).
    ValueError
        ``top`` is not a positive integer or ``None``, or would hide a taxon
        named ``"other"``.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/plotting`

    The bars are ``bt.fn.contributions``' table, drawn as stored: for a
    pathway they need not reach its community value. Plot the stratified
    modality of ``bt.fn.renorm(mdata, "relab")`` for shares of each sample's
    community total.

    Examples
    --------
    >>> import biotapy as bt
    >>> ax = bt.pl.contributions(bt.datasets.toy_humann()["function_by_taxon"], "2.7.1.2")
    >>> len(ax.patches)  # 6 samples x 2 taxa
    12
    """
    strata = function_contributions(adata, function, top=top)
    # fn.contributions adds a column ("other") only when top leaves taxa out; it is drawn grey, as a missing group.
    other = top is not None and strata.shape[1] > top
    n_taxa = strata.shape[1] - 1 if other else strata.shape[1]
    ax = new_axes(ax)
    _stack(ax, strata.to_numpy(), strata.columns.tolist(), colors=_colors(n_taxa, missing=other))
    label_ticks(ax, strata.index.tolist(), axis="x")
    ax.set_xlabel("sample")
    ax.set_ylabel("abundance")
    ax.set_title(function)
    ax.legend(title="taxon")
    return ax


def heatmap(adata: AnnData, *, layer: str | None = None, ax: "Axes | None" = None) -> "Axes":
    """Heatmap of the table, samples as columns and features as rows, on a log colour scale.

    Parameters
    ----------
    adata
        Samples x features, drawn in ``obs`` and ``var`` order.
    layer
        Plot ``layers[layer]``, such as ``"relative"`` from :func:`biotapy.pp.relative`, instead of ``X``.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        The image, features x samples, with a colour bar. Zeros are drawn black.

    Raises
    ------
    KeyError
        ``layers[layer]`` is missing.
    ValueError
        The table holds a negative value, as ``layers["clr"]`` does, or no positive
        value, so the log scale has nothing to show.

    Notes
    -----
    R equivalent: ``phyloseq::plot_heatmap``
    Guide: :doc:`/guide/plotting`

    phyloseq orders samples and features by the angle of their scores on a
    two-axis NMDS that it computes on the fly; biotapy keeps ``obs`` and ``var``
    order and computes nothing, so sort first (see the guide). Colours are
    phyloseq's: ``#000033`` to ``#66CCFF`` on a log scale, zeros black. Tick labels
    are drawn up to 250 names per axis, as phyloseq's ``max.label``. matplotlib needs
    a dense array: the table is densified once (8 bytes x samples x features).

    Examples
    --------
    >>> import biotapy as bt
    >>> ax = bt.pl.heatmap(bt.datasets.toy())
    >>> ax.images[0].get_array().shape  # features x samples
    (8, 6)
    """
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    # matplotlib's imshow needs a dense array (rules.md R6.2): one dense copy, features x samples.
    dense = table(adata, layer, func="pl.heatmap").T.toarray()
    if not (dense > 0).any():
        table_name = "adata: X" if layer is None else f"layer={layer!r}: layers[{layer!r}]"
        msg = f"{table_name} holds no positive value to draw on a log scale"
        raise ValueError(msg)
    ax = new_axes(ax)
    cmap = LinearSegmentedColormap.from_list("phyloseq", [_LOW, _HIGH]).with_extremes(bad=_ZERO)
    image = ax.imshow(dense, cmap=cmap, norm=LogNorm(), aspect="auto", interpolation="nearest")
    ax.figure.colorbar(image, ax=ax, label=layer or "abundance")
    label_ticks(ax, adata.obs_names.tolist(), axis="x")
    label_ticks(ax, adata.var_names.tolist(), axis="y")
    ax.set_xlabel("sample")
    ax.set_ylabel("feature")
    return ax


def _stack(ax: "Axes", heights: npt.NDArray[np.float64], labels: list[str], *, colors: list[RGBA]) -> None:
    """One bar per row of ``heights``, stacking one coloured segment per column in column order."""
    positions, bottom = np.arange(heights.shape[0]), np.zeros(heights.shape[0])
    for column, (label, color) in enumerate(zip(labels, colors, strict=True)):
        ax.bar(positions, heights[:, column], bottom=bottom, color=color, label=label)
        bottom += heights[:, column]


def _segments(
    adata: AnnData, values: sp.csr_matrix, fill: str
) -> tuple[npt.NDArray[np.float64], list[str], list[RGBA]]:
    in_var, in_obs = fill in adata.var.columns, fill in adata.obs.columns
    if in_var == in_obs:
        msg = f"fill={fill!r} must be a column of exactly one of var and obs"
        raise KeyError(msg)
    if in_var:
        codes, labels, colors = groups(cast("pd.Series", adata.var[fill]), arg="fill")
        # Samples x groups: small, unlike X, so it is safe to densify.
        return np.asarray(sum_by(values, codes, len(labels)).toarray(), dtype=np.float64), labels, colors
    codes, labels, colors = obs_groups(adata, fill, arg="fill")
    segments = np.zeros((adata.n_obs, len(labels)))
    segments[np.arange(adata.n_obs), codes] = np.asarray(values.sum(axis=1)).ravel()
    return segments, labels, colors


def _by_x(
    adata: AnnData, segments: npt.NDArray[np.float64], x: str | None
) -> tuple[npt.NDArray[np.float64], list[str]]:
    if x is None:
        return segments, adata.obs_names.tolist()
    codes, names, _ = obs_groups(adata, x, arg="x")
    indicator = np.zeros((len(names), adata.n_obs))
    indicator[codes, np.arange(adata.n_obs)] = 1
    return indicator @ segments, names
