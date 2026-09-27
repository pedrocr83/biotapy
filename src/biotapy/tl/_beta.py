"""Beta diversity: sample x sample distance matrices in ``obsp``."""

from typing import Literal, get_args

import numpy as np
import pandas as pd
from anndata import AnnData
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity

from biotapy._core import TreeData, as_csr, get_skbio_tree, require_counts

BetaMetric = Literal["braycurtis", "jaccard"]

# The call that writes each obsp key, for the error raised when a key is missing.
_WRITTEN_BY = {
    "braycurtis": "bt.tl.beta(adata, metric='braycurtis', inplace=True)",
    "jaccard": "bt.tl.beta(adata, metric='jaccard', inplace=True)",
    "unweighted_unifrac": "bt.tl.unifrac(tdata, inplace=True)",
    "weighted_unifrac": "bt.tl.unifrac(tdata, weighted=True, inplace=True)",
}


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
    features) and the result takes 8 bytes x samples x samples. scikit-bio's
    working copies add to that: measured with tracemalloc on scikit-bio 0.7.4, peak
    memory is the dense copy plus 1.5 results (its condensed and square matrices),
    or two results while the matrix becomes a DataFrame, whichever is larger, and
    ``"jaccard"`` adds a 1-byte presence/absence copy of ``X``.

    References
    ----------
    Bray JR, Curtis JT (1957) An ordination of the upland forest communities of southern
    Wisconsin. Ecological Monographs 27:325-349.

    Jaccard P (1912) The distribution of the flora in the alpine zone. New Phytologist
    11:37-50.

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
    ValueError
        ``weighted=True`` and ``X`` is not raw counts (``x_kind`` ``"counts"`` and whole numbers).

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

    Weighted UniFrac needs raw counts: scikit-bio's tree code casts abundances to
    integers, which truncates proportions to 0. Unweighted UniFrac uses presence
    only and runs on any abundance.

    References
    ----------
    Lozupone C, Knight R (2005) UniFrac: a new phylogenetic method for comparing microbial
    communities. Applied and Environmental Microbiology 71:8228-8235.

    Lozupone CA et al. (2007) Quantitative and qualitative beta diversity measures lead to
    different insights into factors that structure microbial communities. Applied and
    Environmental Microbiology 73:1576-1585.

    Examples
    --------
    >>> import biotapy as bt
    >>> round(float(bt.tl.unifrac(bt.datasets.toy()).loc["s1", "s4"]), 3)
    0.092
    """
    tree = get_skbio_tree(tdata)
    if weighted:
        require_counts(tdata, func="tl.unifrac(weighted=True)")
    # scikit-bio rejects sparse input (rules.md R6.2): one dense copy of X.
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


def stored_distances(adata: AnnData, key: str) -> DistanceMatrix:
    """``obsp[key]`` as a scikit-bio DistanceMatrix; a missing key names the call that writes it."""
    if key not in adata.obsp:
        call = _WRITTEN_BY.get(key, "bt.tl.beta or bt.tl.unifrac with inplace=True")
        msg = f"distance={key!r}: no obsp[{key!r}]; run {call} first"
        raise KeyError(msg)
    values = np.asarray(adata.obsp[key], dtype=np.float64)
    if np.isnan(values).any():
        msg = f"distance={key!r}: obsp[{key!r}] holds NaN, as between two all-zero samples; drop them first"
        raise ValueError(msg)
    return DistanceMatrix(values, ids=adata.obs_names.tolist())
