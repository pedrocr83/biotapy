"""Per-sample transforms: add one layer, keep everything else."""

import numpy as np
from anndata import AnnData

from biotapy._core import add_provenance, as_csr, divide_rows


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
    sums = np.asarray(X.sum(axis=1), dtype=np.float64).ravel()
    out = adata.copy()
    out.layers["relative"] = divide_rows(X, sums)
    add_provenance(out, "pp.relative")
    return out
