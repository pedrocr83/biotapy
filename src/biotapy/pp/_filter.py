"""Filters that keep a subset of features or samples."""

import numpy as np
from anndata import AnnData

from biotapy._core import add_provenance, as_csr, feature_subset


def filter_features(adata: AnnData, *, min_prevalence: float | None = None, min_total: float | None = None) -> AnnData:
    """Keep features that are present in enough samples and have enough reads.

    Parameters
    ----------
    adata
        Samples x features.
    min_prevalence
        Keep features that are non-zero in at least this fraction of samples, from 0 to 1.
    min_total
        Keep features whose total over all samples is at least this.

    Returns
    -------
    AnnData
        Same type as ``adata`` with the kept features in their original order; a
        TreeData keeps their subtree. Both thresholds are inclusive and, when both
        are given, a feature must pass both. ``layers``, ``obsm``, ``obsp``,
        ``varm``, ``varp`` and every non-``biotapy`` ``uns`` key are dropped
        because they described the old features.

    Raises
    ------
    ValueError
        Neither threshold is given, ``min_prevalence`` is outside 0 to 1, or no
        feature passes.

    Notes
    -----
    R equivalent: ``phyloseq::filter_taxa``
    Guide: :doc:`/guide/filtering`

    In R: ``filter_taxa(physeq, function(x) sum(x > 0) >= p * length(x), prune = TRUE)``
    for ``min_prevalence=p``, and ``function(x) sum(x) >= n`` for ``min_total=n``.
    biotapy keeps a feature when ``present / n_obs >= min_prevalence``, while phyloseq's
    ``sum(x > 0) >= p * length(x)`` can drop it at an exact boundary through floating
    point: 7 of 25 samples at ``p = 0.28``, since ``0.28 * 25`` is ``7.000000000000001``.

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.pp.filter_features(bt.datasets.toy(), min_prevalence=1.0).n_vars
    2
    """
    if min_prevalence is None and min_total is None:
        msg = "pass min_prevalence=, min_total= or both"
        raise ValueError(msg)
    X = as_csr(adata.X)
    keep = np.ones(adata.n_vars, dtype=bool)
    if min_prevalence is not None:
        if not 0 <= min_prevalence <= 1:
            msg = f"min_prevalence must be between 0 and 1, got {min_prevalence}"
            raise ValueError(msg)
        # Count stored non-zero values per column: a CSR matrix may also store explicit zeros.
        present = np.bincount(X.indices[X.data != 0], minlength=adata.n_vars)
        # Divide rather than multiply: 7 / 25 >= 0.28 holds, 7 >= 0.28 * 25 does not.
        keep &= present / adata.n_obs >= min_prevalence
    if min_total is not None:
        keep &= np.asarray(X.sum(axis=0)).ravel() >= min_total
    if not keep.any():
        msg = f"no feature passes min_prevalence={min_prevalence}, min_total={min_total}"
        raise ValueError(msg)
    out = feature_subset(adata, np.flatnonzero(keep))
    add_provenance(out, "pp.filter_features", min_prevalence=min_prevalence, min_total=min_total)
    return out


def filter_samples(adata: AnnData, min_depth: float) -> AnnData:
    """Keep samples with at least ``min_depth`` reads.

    Parameters
    ----------
    adata
        Samples x features.
    min_depth
        Keep samples whose total over all features is at least this.

    Returns
    -------
    AnnData
        Same type as ``adata`` with the kept samples in their original order.
        Every slot is kept and subset by AnnData indexing, so ``obsp`` distances
        stay valid; features that are now all-zero are kept. Kept ordinations
        (``obsm`` and the ``pcoa``/``nmds`` summaries in ``uns['biotapy']``) were
        computed with the dropped samples included, so recompute them.

    Raises
    ------
    ValueError
        No sample has ``min_depth`` reads.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/filtering`

    In R: ``prune_samples(sample_sums(physeq) >= min_depth, physeq)``.

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.pp.filter_samples(bt.datasets.toy(), 70).n_obs
    4
    """
    depth = np.asarray(as_csr(adata.X).sum(axis=1)).ravel()
    keep = np.flatnonzero(depth >= min_depth)
    if keep.size == 0:
        msg = f"no sample has min_depth={min_depth} reads; the deepest has {depth.max():g}"
        raise ValueError(msg)
    out = adata[keep].copy()
    add_provenance(out, "pp.filter_samples", min_depth=min_depth)
    return out
