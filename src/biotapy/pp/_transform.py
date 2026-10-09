"""Per-sample transforms: add one layer, keep everything else."""

import numpy as np
from anndata import AnnData
from skbio.stats.composition import clr as skbio_clr

from biotapy._core import add_provenance, as_csr, divide_rows, pseudocounted


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
    8 bytes x samples x features. Peak memory is three to five such arrays (the
    dense copy of ``X``, the transform's temporaries and its output; 3.1x to 4.8x
    measured on 400 x 512 and 2,000 x 2,000 tables, the most when ``X`` is mostly
    non-zero), so budget for that on large tables.

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
    values = pseudocounted(adata.X, pseudocount, func="pp.clr")
    out = adata.copy()
    out.layers["clr"] = skbio_clr(values)
    add_provenance(out, "pp.clr", pseudocount=pseudocount)
    return out
