"""Functional redundancy of each sample from taxon abundances and genome contents (Tian et al. 2020)."""

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.spatial.distance import pdist, squareform

from biotapy._core import as_csr, divide_rows, warn_user

COLUMNS = ["taxonomic_diversity", "functional_diversity", "redundancy", "normalized_redundancy"]


def functional_redundancy(adata: AnnData, *, traits: pd.DataFrame) -> pd.DataFrame:
    r"""Taxonomic diversity, functional diversity and functional redundancy of every sample.

    Parameters
    ----------
    adata
        Samples x taxa abundances whose ``var_names`` are ``traits``' row
        ids, such as an ASV table. Tian et al. use relative organism
        abundances: divide 16S read counts by each ASV's predicted 16S copy
        number first (see the guide).
    traits
        Taxa x genes copy numbers, the genome content of each taxon, such as
        ``bt.io.read_picrust2_traits("EC_predicted.tsv.gz")``. Non-negative
        and finite; row ids unique.

    Returns
    -------
    pandas.DataFrame
        One row per sample (index ``obs_names``) and the columns
        ``taxonomic_diversity`` (Gini-Simpson), ``functional_diversity``
        (Rao's quadratic entropy), ``redundancy`` (their difference) and
        ``normalized_redundancy`` (``redundancy / taxonomic_diversity``). A
        sample with no abundance on a taxon in ``traits`` gets NaN in every
        column; a sample with one such taxon gets 0, 0, 0 and NaN.

    Raises
    ------
    TypeError
        ``adata`` is not an AnnData (for a MuData, pass one modality);
        ``traits`` is not a DataFrame or has a non-numeric column.
    ValueError
        ``traits`` repeats a row id or holds a missing, negative or
        infinite value; ``adata`` repeats a taxon id, its ``X`` holds a
        missing, negative or infinite value, or a sample's total abundance
        overflows; no taxon of ``adata`` has a row in ``traits`` (the
        message shows ids from both, usually an id mismatch).

    Warns
    -----
    UserWarning
        Some taxa of ``adata`` have no row in ``traits`` (PICRUSt2 drops
        ASVs above its NSTI cutoff, and its ``RARE`` group has no genome).
        They are left out, abundance included, and the warning names up to
        three and the count.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/function`

    With :math:`p_i` the relative abundance of taxon :math:`i` over the taxa
    in ``traits``, and :math:`d_{ij}` the weighted Jaccard distance between
    the gene copy numbers of taxa :math:`i` and :math:`j` (Tian et al.,
    Eq. 7), ``taxonomic_diversity`` is :math:`\sum_{i \ne j} p_i p_j`
    (Eq. 2), ``functional_diversity`` is
    :math:`\sum_{i \ne j} d_{ij} p_i p_j` (Eq. 3) and ``redundancy`` is
    :math:`\sum_{i \ne j} (1 - d_{ij}) p_i p_j` (Eqs. 1 and 4). Two taxa
    with no gene at all share nothing, so their distance is 1.

    The distances form one dense taxa x taxa ``float64`` matrix over the
    taxa with abundance in some sample: 8 bytes x taxa², 32 MB for 2,000
    taxa and 800 MB for 10,000. While it is built the peak is about 1.5 x
    8 bytes x taxa² (measured at 3,000 taxa), for the square and the
    condensed half. ``traits`` is read once into a dense ``float64`` copy,
    all its rows, not only the taxa of ``adata``. Time grows as taxa² x
    genes: with 2,500 genes and 100 samples, 2,000 taxa took 5 s and 10,000
    taxa 6 minutes and 1.7 GB beyond the inputs (measured on one machine). On a large ASV table,
    filter rare taxa first with ``bt.pp.filter_features``.

    References
    ----------
    Tian L et al. (2020) Deciphering functional redundancy in the human microbiome.
    Nature Communications 11:6217.

    Examples
    --------
    >>> import anndata as ad
    >>> import numpy as np
    >>> import pandas as pd
    >>> import biotapy as bt
    >>> traits = pd.DataFrame([[2, 1, 0], [1, 1, 1], [0, 0, 3]], index=["A", "B", "C"], columns=["g1", "g2", "g3"])
    >>> adata = ad.AnnData(
    ...     np.array([[2.0, 1.0, 1.0]]), obs=pd.DataFrame(index=["s1"]), var=pd.DataFrame(index=["A", "B", "C"])
    ... )
    >>> bt.fn.functional_redundancy(adata, traits=traits).round(3).loc["s1"].tolist()
    [0.625, 0.475, 0.15, 0.24]
    """
    genomes = _genomes(traits)
    X, taxa = _abundances(adata, genomes.index)
    totals = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    shares = divide_rows(X, totals)
    distance = _weighted_jaccard(genomes.loc[taxa].to_numpy())
    diversity = _quadratic(shares, distance)
    # The similarity 1 - d, built in place over the distances, with i = j left out as in Eq. 4.
    np.subtract(1.0, distance, out=distance)
    np.fill_diagonal(distance, 0.0)
    redundancy = _quadratic(shares, distance)
    # TD as FD + FR (equal to 1 - sum(p**2) up to rounding), so 0 <= FD, FR <= TD hold exactly.
    total = diversity + redundancy
    normalized = np.divide(redundancy, total, out=np.full_like(total, np.nan), where=total > 0)
    out = pd.DataFrame(
        dict(zip(COLUMNS, (total, diversity, redundancy, normalized), strict=True)), index=adata.obs_names
    )
    out.loc[totals == 0] = np.nan
    return out


def _genomes(traits: pd.DataFrame) -> pd.DataFrame:
    """``traits`` as float64 with text row ids, after checking it is a copy-number table."""
    if not isinstance(traits, pd.DataFrame):
        msg = f"traits must be a pandas.DataFrame of taxa x genes, such as bt.io.read_picrust2_traits returns, not {type(traits).__name__}"
        raise TypeError(msg)
    text = [str(column) for column, dtype in traits.dtypes.items() if not pd.api.types.is_numeric_dtype(dtype)]
    if text:
        msg = f"traits must hold copy numbers; non-numeric columns: {text[:3]}"
        raise TypeError(msg)
    genomes = pd.DataFrame(
        traits.to_numpy(dtype=np.float64, na_value=np.nan), index=traits.index.astype(str), columns=traits.columns
    )
    repeated = genomes.index[genomes.index.duplicated()].unique().tolist()
    if repeated:
        msg = f"traits repeats row ids: {repeated[:3]}"
        raise ValueError(msg)
    values = genomes.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any():
        msg = "traits holds a missing, negative or infinite copy number"
        raise ValueError(msg)
    return genomes


def _abundances(adata: AnnData, known: pd.Index) -> tuple[sp.csr_matrix, pd.Index]:
    """``X`` over the taxa that have traits and abundance somewhere, and those taxa; warns about the others."""
    if not isinstance(adata, AnnData):
        msg = (
            f"adata must be an AnnData of samples x taxa, not {type(adata).__name__}; "
            "for a MuData, pass one modality"
        )
        raise TypeError(msg)
    repeated = adata.var_names[adata.var_names.duplicated()].unique().tolist()
    if repeated:
        msg = f"adata repeats taxon ids: {repeated[:3]}"
        raise ValueError(msg)
    X = as_csr(adata.X)
    if not np.isfinite(X.data).all() or (X.data < 0).any():
        msg = "adata: X holds a missing, negative or infinite abundance"
        raise ValueError(msg)
    found = adata.var_names.isin(known)
    if not found.any():
        msg = (
            "no taxon of adata has a row in traits; "
            f"adata: {adata.var_names[:3].tolist()}, traits: {known[:3].tolist()}"
        )
        raise ValueError(msg)
    if not found.all():
        missing = adata.var_names[~found]
        warn_user(
            f"fn.functional_redundancy: {len(missing)} of {adata.n_vars} taxa have no row in traits "
            f"and are left out, with their abundance: {missing[:3].tolist()}"
        )
    keep = np.flatnonzero(found & (np.asarray(X.sum(axis=0)).ravel() > 0))
    X = X[:, keep]
    with np.errstate(over="ignore"):
        overflowing = ~np.isfinite(np.asarray(X.sum(axis=1, dtype=np.float64)).ravel())
    if overflowing.any():
        msg = f"adata: the total abundance of a sample overflows: {adata.obs_names[overflowing][:3].tolist()}"
        raise ValueError(msg)
    return X, adata.var_names[keep]


def _weighted_jaccard(genomes: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Tian et al.'s Eq. 7 for every pair of rows, as a square matrix with a zero diagonal."""
    if genomes.shape[0] < 2:
        return np.zeros((genomes.shape[0], genomes.shape[0]))
    # For non-negative vectors, 1 - sum(min)/sum(max) = 2 BC / (1 + BC), with BC the Bray-Curtis dissimilarity.
    bray_curtis = pdist(genomes, "braycurtis")
    # In place on the condensed array, so only one other condensed-sized temporary exists at a time.
    np.divide(bray_curtis, 1 + bray_curtis, out=bray_curtis)
    bray_curtis *= 2
    # Two taxa without any gene give 0/0 (SciPy returns NaN): they share nothing, so their distance is 1.
    bray_curtis[np.isnan(bray_curtis)] = 1.0
    return squareform(bray_curtis)


def _quadratic(shares: sp.csr_matrix, matrix: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """``p @ matrix @ p`` for every sample's row ``p`` of ``shares``."""
    return np.asarray(shares.multiply(shares @ matrix).sum(axis=1), dtype=np.float64).ravel()
