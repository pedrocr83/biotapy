"""Aggregation along the taxonomy."""

from typing import cast

import numpy as np
import pandas as pd
from anndata import AnnData

from biotapy._core import add_provenance, argmax_by, as_csr, feature_subset, split_ranks, sum_by


def tax_glom(adata: AnnData, rank: str, *, dropna: bool = True) -> AnnData:
    """Aggregate features to a taxonomic rank.

    Features that share a lineage down to ``rank`` are summed into one
    feature, represented by the group's most abundant member.

    Parameters
    ----------
    adata
        Samples x features with taxonomy columns in ``var``.
    rank
        Canonical rank name, for example ``"genus"``.
    dropna
        Drop features with no value at ``rank`` before aggregating.

    Returns
    -------
    AnnData
        Same type as ``adata``; a TreeData keeps the representatives' subtree.
        Ranks below ``rank`` are ``NaN``. ``layers``, ``obsm`` and ``obsp`` are
        dropped because they described the old features.

    Raises
    ------
    KeyError
        ``rank`` is not a taxonomy column.
    ValueError
        No feature has a value at ``rank``.

    Notes
    -----
    R equivalent: ``phyloseq::tax_glom``, ``mia::agglomerateByRank``
    Guide: :doc:`/guide/aggregation`

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.pp.tax_glom(bt.datasets.toy(), "phylum").n_vars
    3
    """
    upto, below = split_ranks(adata, rank)
    # anndata types .var as DataFrame | Dataset2D (its lazy/backed variant); biotapy's
    # data model (data-model-slots) guarantees a real DataFrame here.
    var_upto = cast("pd.DataFrame", adata.var[upto])
    # phyloseq keys groups by paste(lineage, collapse=";_;"), NA included
    lineage = var_upto.astype("string").fillna("NA").agg(";_;".join, axis=1)
    if dropna:
        lineage = lineage.where(adata.var[rank].notna())
    codes, groups = pd.factorize(lineage)
    if groups.size == 0:
        msg = f"no feature has a value at rank={rank!r}"
        raise ValueError(msg)
    X = as_csr(adata.X)
    keep = argmax_by(np.asarray(X.sum(axis=0), dtype=np.float64).ravel(), codes)
    order = np.argsort(keep)
    out = feature_subset(adata, keep[order])
    out.X = sum_by(X, codes, groups.size)[:, order]
    out.var[below] = np.nan
    add_provenance(out, "pp.tax_glom", rank=rank, dropna=dropna)
    return out
