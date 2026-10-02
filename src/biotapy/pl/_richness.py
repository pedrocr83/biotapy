"""Alpha diversity per sample, as stored by ``tl.alpha``."""

from typing import TYPE_CHECKING

import numpy as np
from anndata import AnnData

from ._common import label_ticks, new_axes, obs_groups, scatter

if TYPE_CHECKING:
    from matplotlib.axes import Axes


def richness(
    adata: AnnData, metric: str, *, x: str | None = None, color: str | None = None, ax: "Axes | None" = None
) -> "Axes":
    """One point per sample for a stored alpha diversity metric.

    Parameters
    ----------
    adata
        Samples x features with ``obs['alpha_<metric>']``, written by
        :func:`biotapy.tl.alpha` with ``inplace=True``.
    metric
        The metric to plot, such as ``"shannon"``.
    x
        An ``obs`` column whose values place the points. By default one position per sample.
    color
        An ``obs`` column that colours the points, with a legend.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        The points, one per sample with a finite value.

    Raises
    ------
    KeyError
        ``obs['alpha_<metric>']`` is missing (the message names the call that writes
        it), or ``x`` or ``color`` is not an ``obs`` column.
    TypeError
        ``x`` or ``color`` is a numeric column.

    Notes
    -----
    R equivalent: ``phyloseq::plot_richness``
    Guide: :doc:`/guide/plotting`

    phyloseq computes ``estimate_richness`` and draws every measure in its own
    facet; biotapy reads what :func:`biotapy.tl.alpha` stored and draws one metric
    per axes. Samples with a NaN value, such as Shannon of an all-zero sample, are
    left out, as ``geom_point(na.rm = TRUE)`` does. There are no standard-error bars:
    ``tl.alpha`` stores none.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)
    >>> ax = bt.pl.richness(tdata, "shannon", x="group")
    >>> ax.collections[0].get_offsets().shape
    (6, 2)
    """
    column = f"alpha_{metric}"
    if column not in adata.obs.columns:
        msg = f"metric={metric!r}: no obs[{column!r}]; run bt.tl.alpha(adata, metrics=[{metric!r}], inplace=True) first"
        raise KeyError(msg)
    values = np.asarray(adata.obs[column], dtype=np.float64)
    by = obs_groups(adata, color, arg="color") if color is not None else None
    if x is None:
        positions, names = np.arange(adata.n_obs, dtype=np.float64), adata.obs_names.tolist()
    else:
        codes, names, _ = obs_groups(adata, x, arg="x")
        positions = codes.astype(np.float64)
    finite = np.isfinite(values)
    if by is not None:
        by = (by[0][finite], by[1], by[2])
    ax = new_axes(ax)
    scatter(ax, positions[finite], values[finite], by=by, title=color)
    label_ticks(ax, names, axis="x")
    ax.set_xlabel(x or "sample")
    ax.set_ylabel(metric)
    return ax
