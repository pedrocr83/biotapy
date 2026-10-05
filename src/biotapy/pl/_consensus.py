"""Where differential abundance methods agree: a dot per feature and method, from da.consensus's table."""

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ._common import MISSING_COLOR, label_ticks, new_axes

if TYPE_CHECKING:
    from matplotlib.axes import Axes

# tab10's red and blue: a call's sign. A hollow grey dot is a feature the method tested but did not call.
_KINDS = {
    "effect > 0": {"color": "#d62728"},
    "effect < 0": {"color": "#1f77b4"},
    "not significant": {"facecolors": "none", "edgecolors": MISSING_COLOR},
}


def _check(table: pd.DataFrame, *, top: int) -> list[str]:
    """The method names in ``table``, after checking it is a consensus table and ``top`` a positive int."""
    if not isinstance(table, pd.DataFrame):
        msg = f"table must be a pandas DataFrame, got {type(table).__name__}"
        raise TypeError(msg)
    methods = [column.removeprefix("significant_") for column in table.columns if column.startswith("significant_")]
    needed = ["n_significant", "consensus", *[f"{kind}_{m}" for m in methods for kind in ("effect", "qvalue")]]
    absent = [column for column in needed if column not in table.columns]
    if not methods or absent:
        msg = f"table lacks {absent or ['significant_<method>']}; pl.consensus draws the table bt.da.consensus returns"
        raise KeyError(msg)
    if isinstance(top, bool) or not isinstance(top, int | np.integer):
        msg = f"top must be an integer, got {top!r}"
        raise TypeError(msg)
    if top < 1:
        msg = f"top must be at least 1, got {top}"
        raise ValueError(msg)
    return methods


def consensus(table: pd.DataFrame, *, top: int = 30, ax: "Axes | None" = None) -> "Axes":
    """A dot matrix of which methods call which features, from :func:`biotapy.da.consensus`.

    Parameters
    ----------
    table
        The table :func:`biotapy.da.consensus` returns.
    top
        Draw at most this many features: those called by the most methods, then with
        the largest mean absolute effect, then in table order. A ``top`` above 30
        needs a taller figure, passed as ``ax``.
    ax
        Axes to draw on; by default a new figure's.

    Returns
    -------
    matplotlib.axes.Axes
        One row per feature (top row first) and one column per method: a filled dot
        coloured by the sign of the effect where the method calls the feature, a
        hollow dot where it tested the feature without calling it, nothing where it did
        not test it. Consensus features have bold labels.

    Raises
    ------
    TypeError
        ``table`` is not a DataFrame, or ``top`` is not an integer.
    KeyError
        ``table`` lacks the columns :func:`biotapy.da.consensus` writes.
    ValueError
        ``top`` is below 1, or no method calls any feature.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/differential_abundance`

    One axes, as every ``pl`` function returns; an UpSet plot of the same calls
    needs two panels. For a two-level ``group`` the effects of all methods are log2
    fold changes, which makes their mean comparable when ranking features. For a
    numeric ``group`` they are not (``da.ancombc2`` gives the change per unit,
    ``da.linda`` per standard deviation), so the ranking mixes units: rank by
    ``n_significant`` and read each method's effect on its own.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> table = bt.da.consensus([bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")])
    >>> ax = bt.pl.consensus(table)
    >>> [label.get_text() for label in ax.get_yticklabels()]
    ['f6', 'f7', 'f8']
    """
    methods = _check(table, top=top)
    called = table[table["n_significant"] > 0]
    if called.empty:
        msg = "no method calls any feature significant: nothing to draw"
        raise ValueError(msg)
    strength = called[[f"effect_{m}" for m in methods]].abs().mean(axis=1).to_numpy()
    rows = called.iloc[np.lexsort((-strength, -called["n_significant"].to_numpy()))[:top]]
    significant = rows[[f"significant_{m}" for m in methods]].to_numpy(bool)
    effect = rows[[f"effect_{m}" for m in methods]].to_numpy(np.float64)
    tested = rows[[f"qvalue_{m}" for m in methods]].notna().to_numpy()
    y, x = np.indices(significant.shape)
    masks = [significant & (effect > 0), significant & (effect < 0), tested & ~significant]
    ax = new_axes(ax)
    for (label, style), mask in zip(_KINDS.items(), masks, strict=True):
        if mask.any():
            ax.scatter(x[mask], y[mask], label=label, **style)  # type: ignore[arg-type]
    ax.set_xticks(np.arange(len(methods)), labels=methods)
    label_ticks(ax, rows.index.astype(str).tolist(), axis="y")
    for tick, bold in zip(ax.get_yticklabels(), rows["consensus"].to_numpy(bool), strict=False):
        tick.set_fontweight("bold" if bold else "normal")
    ax.set_xlim(-0.5, len(methods) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    # One row above the axes: inside they cover dots, beside them a default save clips them.
    kinds = len(ax.get_legend_handles_labels()[0])
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=kinds, frameon=False, borderaxespad=0)
    return ax
