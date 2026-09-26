"""Per-sample transforms: add one layer, keep everything else."""

from typing import cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import add_provenance, as_csr


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
    # anndata's stub types X as _XDataType | None (backed datasets, missing X); biotapy's
    # contract (data-model-slots) guarantees X is populated and array-like here.
    X = as_csr(cast("sp.spmatrix | sp.sparray | npt.ArrayLike", adata.X))
    sums = np.asarray(X.sum(axis=1), dtype=np.float64).ravel()
    scale = np.divide(1.0, sums, out=np.zeros_like(sums), where=sums > 0)
    out = adata.copy()
    out.layers["relative"] = sp.csr_matrix(sp.diags(scale) @ X)
    add_provenance(out, "pp.relative")
    return out
