"""Beta diversity: sample x sample distance matrices in ``obsp``."""

from typing import Literal, get_args

import pandas as pd
from anndata import AnnData
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity

from biotapy._core import TreeData, as_csr, get_skbio_tree

BetaMetric = Literal["braycurtis", "jaccard"]


def beta(adata: AnnData, *, metric: BetaMetric = "braycurtis", inplace: bool = False) -> pd.DataFrame | None:
    """Distances between every pair of samples.

    Parameters
    ----------
    adata
        Samples x features.
    metric
        ``"braycurtis"``, or ``"jaccard"`` on presence/absence.
    inplace
        Write the matrix to ``obsp[metric]`` and return ``None``.

    Returns
    -------
    pandas.DataFrame or None
        Symmetric samples x samples distances with a zero diagonal, indexed by
        ``obs_names``. Two all-zero samples are NaN apart under Bray-Curtis and 0
        apart under Jaccard (scikit-bio's convention; vegan's binary Jaccard gives NaN).

    Raises
    ------
    ValueError
        ``metric`` is not ``"braycurtis"`` or ``"jaccard"``.

    Notes
    -----
    R equivalent: ``phyloseq::distance``
    Guide: :doc:`/guide/diversity`

    ``"jaccard"`` matches ``phyloseq::distance(physeq, "jaccard", binary = TRUE)``:
    without ``binary = TRUE``, vegan computes a quantitative Jaccard instead.
    scikit-bio needs dense input, so ``X`` is densified once (8 bytes x samples x
    features) and the result takes 8 bytes x samples x samples.

    Examples
    --------
    >>> import biotapy as bt
    >>> round(float(bt.tl.beta(bt.datasets.toy()).loc["s1", "s4"]), 3)
    0.708
    """
    if metric not in get_args(BetaMetric):
        msg = f"metric must be one of {list(get_args(BetaMetric))}, got {metric!r}"
        raise ValueError(msg)
    # scikit-bio rejects sparse input (rules.md R6.2): one dense copy of X.
    distances = beta_diversity(metric, as_csr(adata.X).toarray(), ids=adata.obs_names.tolist())
    return _store(adata, distances, key=metric, inplace=inplace)


def unifrac(
    tdata: TreeData, *, weighted: bool = False, normalized: bool = True, inplace: bool = False
) -> pd.DataFrame | None:
    """UniFrac distances between every pair of samples.

    Parameters
    ----------
    tdata
        Samples x features with the phylogeny in ``vart['phylo']``.
    weighted
        Weight branches by abundance (weighted UniFrac) instead of presence.
    normalized
        Scale weighted UniFrac to 0-1, as phyloseq does. Ignored when unweighted.
    inplace
        Write the matrix to ``obsp['unweighted_unifrac']`` or
        ``obsp['weighted_unifrac']`` and return ``None``.

    Returns
    -------
    pandas.DataFrame or None
        Symmetric samples x samples distances with a zero diagonal, indexed by ``obs_names``.
        Two all-zero samples are 0 apart.

    Raises
    ------
    TypeError
        ``tdata`` is an AnnData that is not a TreeData.
    KeyError
        ``tdata`` has no ``vart['phylo']``.

    Notes
    -----
    R equivalent: ``phyloseq::UniFrac``
    Guide: :doc:`/guide/diversity`

    scikit-bio needs a root with at most two children. A root with more keeps its
    first child and gets the others under one new zero-length branch, so the tree
    is used rooted where it is drawn; phyloseq instead roots an unrooted tree at a
    random tip. scikit-bio's weighted UniFrac is unnormalized by default; biotapy
    passes ``normalized`` explicitly and defaults to phyloseq's ``TRUE``.
    ``X`` is densified once (8 bytes x samples x features).

    Examples
    --------
    >>> import biotapy as bt
    >>> round(float(bt.tl.unifrac(bt.datasets.toy()).loc["s1", "s4"]), 3)
    0.092
    """
    tree = get_skbio_tree(tdata)
    counts, ids, taxa = as_csr(tdata.X).toarray(), tdata.obs_names.tolist(), tdata.var_names.tolist()
    if weighted:
        distances = beta_diversity("weighted_unifrac", counts, ids=ids, taxa=taxa, tree=tree, normalized=normalized)
    else:
        distances = beta_diversity("unweighted_unifrac", counts, ids=ids, taxa=taxa, tree=tree)
    return _store(tdata, distances, key="weighted_unifrac" if weighted else "unweighted_unifrac", inplace=inplace)


def _store(adata: AnnData, distances: DistanceMatrix, *, key: str, inplace: bool) -> pd.DataFrame | None:
    frame = distances.to_data_frame()
    if not inplace:
        return frame
    adata.obsp[key] = frame.to_numpy()
    return None
