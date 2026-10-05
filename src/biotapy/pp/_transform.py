"""Per-sample transforms: add one layer, keep everything else."""

import numpy as np
import numpy.typing as npt
from anndata import AnnData
from skbio.stats.composition import clr as skbio_clr

from biotapy._core import add_provenance, as_csr, divide_rows, warn_user


def relative(adata: AnnData) -> AnnData:
    """Add per-sample relative abundance as ``layers['relative']``.

    Parameters
    ----------
    adata
        Samples x features; ``X`` holds counts or another non-negative abundance.

    Returns
    -------
    AnnData
        A copy of ``adata`` (a TreeData stays a TreeData) with
        ``layers['relative']``; ``X`` is unchanged.

    Notes
    -----
    R equivalent: ``phyloseq::transform_sample_counts``, ``mia::transformAssay``
    Guide: :doc:`/guide/transforms`

    All-zero samples stay all-zero, where phyloseq returns ``NaN``.

    Examples
    --------
    >>> import biotapy as bt
    >>> out = bt.pp.relative(bt.datasets.toy())
    >>> round(float(out.layers["relative"][0].sum()), 6)
    1.0
    """
    X = as_csr(adata.X)
    sums = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    out = adata.copy()
    out.layers["relative"] = divide_rows(X, sums)
    add_provenance(out, "pp.relative")
    return out


def clr(adata: AnnData, *, pseudocount: float = 0.5) -> AnnData:
    """Add the centred log-ratio transform of each sample as ``layers['clr']``.

    Parameters
    ----------
    adata
        Samples x features; ``X`` holds counts or another non-negative abundance.
    pseudocount
        Added to every value of ``X`` before the logarithm, so zeros have one.

    Returns
    -------
    AnnData
        A copy of ``adata`` (a TreeData stays a TreeData) with ``layers['clr']``,
        a dense float64 array in which every sample sums to 0; ``X`` is unchanged.

    Raises
    ------
    TypeError
        ``pseudocount`` is a bool or not a real number.
    ValueError
        ``pseudocount`` is negative or not finite; ``X`` holds a negative or
        non-finite value; or ``X`` plus ``pseudocount`` still holds a zero.

    Warns
    -----
    UserWarning
        ``pseudocount`` is larger than the smallest non-zero value in ``X``, as when
        the default 0.5 meets relative abundances.

    Notes
    -----
    R equivalent: ``mia::transformAssay``, ``vegan::decostand``
    Guide: :doc:`/guide/transforms`

    Equals ``mia::transformAssay(tse, method = "clr", pseudocount = 0.5)`` and
    ``vegan::decostand(x, "clr", pseudocount = 0.5)``: the pseudocount is added to
    every value, zeros or not, and the logarithm is natural. mia's default
    ``pseudocount = FALSE`` fails on zeros; biotapy defaults to 0.5, the value
    LinDA uses. The transform is scale invariant, so counts and their relative
    abundances give the same result only when the pseudocount is scaled with them.
    An all-zero sample gives an all-zero CLR row.

    CLR has no zeros, so ``X`` is densified once and ``layers['clr']`` is dense:
    8 bytes x samples x features. Peak memory is about three such arrays (the
    dense copy of ``X``, the transform's temporaries and its output; 3.4x measured
    on a 400 x 500 table), so budget for that on large tables.

    References
    ----------
    Aitchison J (1986) The Statistical Analysis of Compositional Data. Chapman & Hall.

    Examples
    --------
    >>> import biotapy as bt
    >>> out = bt.pp.clr(bt.datasets.toy())
    >>> round(float(out.layers["clr"][0, 0]), 3)
    1.048
    """
    values = pseudocounted(adata, pseudocount, func="pp.clr")
    out = adata.copy()
    out.layers["clr"] = skbio_clr(values)
    add_provenance(out, "pp.clr", pseudocount=pseudocount)
    return out


def pseudocounted(
    adata: AnnData, pseudocount: float, *, func: str, columns: npt.NDArray[np.intp] | None = None
) -> npt.NDArray[np.float64]:
    """``X`` (its ``columns``, in that order, if given) as a dense float64 array plus ``pseudocount``, checked > 0."""
    if isinstance(pseudocount, bool) or not isinstance(pseudocount, int | float | np.integer | np.floating):
        msg = f"pseudocount must be a real number, got {pseudocount!r}"
        raise TypeError(msg)
    if not np.isfinite(pseudocount) or pseudocount < 0:
        msg = f"pseudocount must be a finite number >= 0, got {pseudocount!r}"
        raise ValueError(msg)
    X = as_csr(adata.X).astype(np.float64)
    if columns is not None:
        # Reordered while sparse, so the dense copy below is the only one.
        X = X[:, columns]
    if not np.all(np.isfinite(X.data)) or np.any(X.data < 0):
        msg = f"{func} needs finite, non-negative values in X"
        raise ValueError(msg)
    positive = X.data[X.data > 0]
    if positive.size and pseudocount > positive.min():
        warn_user(
            f"pseudocount={pseudocount} is larger than the smallest non-zero value in X ({positive.min():.3g}), "
            "so it swamps the rarest features; for relative abundances pass a pseudocount on their scale"
        )
    # scikit-bio's log-ratio functions need dense input (rules.md R6.2): one dense copy of X.
    values = X.toarray()
    values += pseudocount
    if np.any(values <= 0):
        msg = f"X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to {func}"
        raise ValueError(msg)
    return values
