"""Abundance plots: stacked bars and a heatmap of the table."""

from typing import TYPE_CHECKING, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import sum_by

from ._common import RGBA, groups, label_ticks, new_axes, obs_groups, table

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
    values = table(adata, layer)
    segments, labels, colors = _segments(adata, values, fill)
    heights, names = _by_x(adata, segments, x)
    ax = new_axes(ax)
    positions, bottom = np.arange(len(names)), np.zeros(len(names))
    for column, (label, color) in enumerate(zip(labels, colors, strict=True)):
        ax.bar(positions, heights[:, column], bottom=bottom, color=color, label=label)
        bottom += heights[:, column]
    label_ticks(ax, names, axis="x")
    ax.set_xlabel(x or "sample")
    ax.set_ylabel(layer or "abundance")
    ax.legend(title=fill)
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
        The table holds no positive value, so the log scale has nothing to show.

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
    dense = table(adata, layer).T.toarray()
    if not (dense > 0).any():
        msg = "the table holds no positive value to draw on a log scale"
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
