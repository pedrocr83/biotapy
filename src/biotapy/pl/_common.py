"""Helpers shared by the pl topic files: new axes, groups of a column, colours, the plotted table."""

from typing import TYPE_CHECKING, Literal, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import as_csr, require_categorical

if TYPE_CHECKING:
    from matplotlib.axes import Axes

RGBA = tuple[float, float, float, float]
# One integer code per sample or feature, the group labels, and one colour per label.
Groups = tuple[npt.NDArray[np.intp], list[str], list[RGBA]]

# ggplot2's defaults, which phyloseq's plots inherit: missing values are their own group, "NA", drawn grey50 and last.
MISSING_LABEL = "NA"
MISSING_COLOR = "#7f7f7f"
# phyloseq's plot_heatmap(max.label = 250): past this many names, tick labels overlap into a solid block.
MAX_LABELS = 250
# The call that writes each layer pl can read, for the error raised when it is missing.
_LAYER_WRITTEN_BY = {"relative": "adata = bt.pp.relative(adata)"}


def new_axes(ax: "Axes | None") -> "Axes":
    """``ax``, or the axes of a new pyplot figure, which a notebook displays and ``plt.show()`` shows."""
    if ax is not None:
        return ax
    # Imported here: pyplot costs about 0.4 s, paid by the first plot and never by `import biotapy`.
    import matplotlib.pyplot as plt

    return plt.figure().add_subplot()


def table(adata: AnnData, layer: str | None) -> sp.csr_matrix:
    """``X``, or ``layers[layer]``; a missing layer names the call that writes it."""
    if layer is None:
        return as_csr(adata.X)
    if layer not in adata.layers:
        call = _LAYER_WRITTEN_BY.get(layer, "the function that writes it")
        msg = f"layer={layer!r}: no layers[{layer!r}]; run {call} first"
        raise KeyError(msg)
    return as_csr(adata.layers[layer])


def obs_groups(adata: AnnData, column: str, *, arg: str) -> Groups:
    """Groups of ``obs[column]``; ``arg`` names the parameter in errors."""
    if column not in adata.obs.columns:
        msg = f"{arg}={column!r} is not a column of obs"
        raise KeyError(msg)
    # anndata types obs columns as Series | DataArray (its lazy variant); the data model guarantees a Series.
    return groups(cast("pd.Series", adata.obs[column]), arg=arg)


def groups(values: pd.Series, *, arg: str) -> Groups:
    """Codes, labels and colours: a categorical's order, else sorted values, with missing values last as NA."""
    require_categorical(values, arg=arg, purpose="plots group by category")
    categorical = pd.Categorical(values).remove_unused_categories()
    codes = categorical.codes.astype(np.intp)
    labels = [str(category) for category in categorical.categories]
    missing = bool((codes < 0).any())
    if missing:
        codes = np.where(codes < 0, len(labels), codes)
        labels.append(MISSING_LABEL)
    return codes, labels, _colors(len(labels) - missing, missing=missing)


def _colors(n: int, *, missing: bool) -> list[RGBA]:
    from matplotlib import colormaps
    from matplotlib.colors import to_rgba

    # tab10 and tab20 as long as they last; past 20 groups, n evenly spaced colours, as ggplot2's hue scale.
    cmap = colormaps["tab10"] if n <= 10 else colormaps["tab20"] if n <= 20 else colormaps["turbo"].resampled(n)
    colors = [cmap(i) for i in range(n)]
    return [*colors, to_rgba(MISSING_COLOR)] if missing else colors


def scatter(
    ax: "Axes", x: npt.NDArray[np.float64], y: npt.NDArray[np.float64], *, by: Groups | None, title: str | None
) -> None:
    """One scatter, or one per group that has points, with a legend titled ``title``."""
    if by is None:
        ax.scatter(x, y)
        return
    codes, labels, colors = by
    for code, (label, color) in enumerate(zip(labels, colors, strict=True)):
        keep = codes == code
        # A group whose points were all left out (NaN) gets no legend entry; the others keep their colours.
        if keep.any():
            ax.scatter(x[keep], y[keep], color=color, label=label)
    ax.legend(title=title)


def label_ticks(ax: "Axes", names: list[str], *, axis: Literal["x", "y"]) -> None:
    """One tick per name, or none past ``MAX_LABELS`` names."""
    shown = names if len(names) <= MAX_LABELS else []
    if axis == "x":
        ax.set_xticks(np.arange(len(shown)), labels=shown, rotation=90)
    else:
        ax.set_yticks(np.arange(len(shown)), labels=shown)
