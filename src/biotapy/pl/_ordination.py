"""Ordination plots: samples on two axes, and the scree of a PCoA."""

from typing import TYPE_CHECKING, Literal, get_args

import numpy as np
import numpy.typing as npt
from anndata import AnnData

from ._common import new_axes, obs_groups, scatter

if TYPE_CHECKING:
    from matplotlib.axes import Axes

Basis = Literal["pcoa", "nmds"]
# The call that writes obsm["X_<basis>"] and uns["biotapy"][<basis>], for the error raised when they are missing.
_WRITTEN_BY = {"pcoa": "bt.tl.pcoa(adata, inplace=True)", "nmds": "bt.tl.nmds(adata, inplace=True)"}


def ordination(
    adata: AnnData,
    *,
    basis: Basis = "pcoa",
    components: tuple[int, int] = (1, 2),
    color: str | None = None,
    ax: "Axes | None" = None,
) -> "Axes":
    """Samples on two axes of a stored ordination.

    Parameters
    ----------
    adata
        Samples x features with ``obsm['X_pcoa']`` or ``obsm['X_nmds']``, written by
        :func:`biotapy.tl.pcoa` or :func:`biotapy.tl.nmds` with ``inplace=True``.
    basis
        ``"pcoa"`` or ``"nmds"``.
    components
        The two axes to draw, counted from 1.
    color
        An ``obs`` column that colours the samples, with a legend.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        One point per sample. PCoA axis labels carry the percentage each axis
        explains; an NMDS plot notes the stress in its lower right corner.

    Raises
    ------
    ValueError
        ``basis`` is not ``"pcoa"`` or ``"nmds"``, or ``components`` does not name two stored axes.
    KeyError
        The ordination is missing (the message names the call that writes it), or
        ``color`` is not an ``obs`` column.
    TypeError
        ``color`` is a numeric column.

    Notes
    -----
    R equivalent: ``phyloseq::plot_ordination``
    Guide: :doc:`/guide/plotting`

    As in phyloseq, PCoA labels append ``round(100 * Relative_eig, 1)`` and NMDS
    labels carry no percentage. The stress note is biotapy's; phyloseq shows none.
    Taxa, biplot and split plots are not in 0.4.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> bt.tl.beta(tdata, inplace=True)
    >>> bt.tl.pcoa(tdata, inplace=True)
    >>> bt.pl.ordination(tdata, color="group").get_xlabel()
    'PC1 [93.8%]'
    """
    coords, summary = _stored(adata, basis)
    first, second = _components(components, coords.shape[1])
    by = obs_groups(adata, color, arg="color") if color is not None else None
    ax = new_axes(ax)
    scatter(ax, coords[:, first], coords[:, second], by=by, title=color)
    if basis == "pcoa":
        shares = np.asarray(summary["proportion_explained"], dtype=np.float64)
        ax.set_xlabel(f"PC{first + 1} [{100 * shares[first]:.1f}%]")
        ax.set_ylabel(f"PC{second + 1} [{100 * shares[second]:.1f}%]")
    else:
        ax.set_xlabel(f"NMDS{first + 1}")
        ax.set_ylabel(f"NMDS{second + 1}")
        stress = float(np.asarray(summary["stress"]))
        ax.text(0.99, 0.01, f"stress {stress:.3f}", transform=ax.transAxes, ha="right", va="bottom")
    return ax


def scree(adata: AnnData, *, ax: "Axes | None" = None) -> "Axes":
    """Bars of the proportion of variation each stored PCoA axis explains.

    Parameters
    ----------
    adata
        Samples x features with ``uns['biotapy']['pcoa']``, written by
        :func:`biotapy.tl.pcoa` with ``inplace=True``.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        One bar per stored axis: ``PC1``, ``PC2``, ...

    Raises
    ------
    KeyError
        The PCoA is missing; the message names the call that writes it.

    Notes
    -----
    R equivalent: ``phyloseq::plot_scree``
    Guide: :doc:`/guide/plotting`

    The bars are ``proportion_explained``, each eigenvalue over their sum, as
    phyloseq's ``x / sum(x)``. phyloseq draws every axis; biotapy draws the axes
    ``tl.pcoa`` stored, 10 by default, so pass ``n_components`` for more.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> bt.tl.beta(tdata, inplace=True)
    >>> bt.tl.pcoa(tdata, inplace=True)
    >>> len(bt.pl.scree(tdata).patches)  # at most n_obs - 1 axes
    5
    """
    _, summary = _stored(adata, "pcoa")
    shares = np.asarray(summary["proportion_explained"], dtype=np.float64)
    ax = new_axes(ax)
    positions = np.arange(shares.size)
    ax.bar(positions, shares)
    ax.set_xticks(positions, labels=[f"PC{i}" for i in range(1, shares.size + 1)])
    ax.set_xlabel("axis")
    ax.set_ylabel("proportion explained")
    return ax


def _stored(adata: AnnData, basis: str) -> tuple[npt.NDArray[np.float64], dict[str, object]]:
    if basis not in get_args(Basis):
        msg = f"basis must be one of {list(get_args(Basis))}, got {basis!r}"
        raise ValueError(msg)
    key = f"X_{basis}"
    summary = adata.uns.get("biotapy", {}).get(basis)
    if key not in adata.obsm or summary is None:
        msg = f"basis={basis!r}: no obsm[{key!r}] or uns['biotapy'][{basis!r}]; run {_WRITTEN_BY[basis]} first"
        raise KeyError(msg)
    return np.asarray(adata.obsm[key], dtype=np.float64), dict(summary)


def _components(components: tuple[int, int], n_axes: int) -> tuple[int, int]:
    if len(components) != 2 or not all(1 <= axis <= n_axes for axis in components):
        msg = f"components must be two axis numbers from 1 to {n_axes}, got {components!r}"
        raise ValueError(msg)
    return components[0] - 1, components[1] - 1
