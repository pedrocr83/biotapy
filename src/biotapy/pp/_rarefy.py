"""Rarefaction: subsample every sample to the same depth."""

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from anndata import AnnData
from skbio.stats import subsample_counts

from biotapy._core import add_provenance, as_csr, as_generator, feature_subset, require_counts, warn_user


def rarefy(adata: AnnData, *, depth: int | None = None, seed: int | np.random.Generator | None = None) -> AnnData:
    """Subsample every sample to ``depth`` reads, without replacement.

    Parameters
    ----------
    adata
        Samples x features with raw counts in ``X``.
    depth
        Reads to keep per sample. Defaults to the smallest non-zero sample depth.
    seed
        Seed or generator for the subsampling.

    Returns
    -------
    AnnData
        Same type as ``adata``. Samples with fewer than ``depth`` reads are dropped,
        with one warning naming them. Every kept sample sums to ``depth``. Features
        left all-zero are dropped, a TreeData keeps the remaining features' subtree,
        and ``layers``, ``obsm``, ``obsp``, ``varm``, ``varp`` and every
        non-``biotapy`` ``uns`` key are dropped.

    Raises
    ------
    TypeError
        ``depth`` is given and is not an integer.
    ValueError
        ``X`` does not hold counts, ``depth`` is below 1, or no sample has ``depth`` reads.

    Warns
    -----
    UserWarning
        Samples with fewer than ``depth`` reads were dropped.

    Notes
    -----
    R equivalent: ``phyloseq::rarefy_even_depth``
    Guide: :doc:`/guide/filtering`

    Draws are without replacement (scikit-bio ``subsample_counts``), so no count
    exceeds the original; phyloseq defaults to ``replace = TRUE``. The default depth
    skips all-zero samples, where phyloseq's ``min(sample_sums(physeq))`` would be 0.
    R and NumPy random generators differ, so the counts never match phyloseq's.

    Examples
    --------
    >>> import biotapy as bt
    >>> out = bt.pp.rarefy(bt.datasets.toy(), depth=60, seed=0)
    >>> out.X.sum(axis=1).A1.tolist()
    [60, 60, 60, 60, 60, 60]
    """
    require_counts(adata, func="pp.rarefy")
    if depth is not None and (isinstance(depth, bool) or not isinstance(depth, int | np.integer)):
        msg = f"depth= must be an integer, got {depth}"
        raise TypeError(msg)
    X = as_csr(adata.X)
    sums = np.asarray(X.sum(axis=1)).ravel()
    depth = _smallest_nonzero(sums) if depth is None else int(depth)
    if depth < 1:
        msg = f"depth must be at least 1, got {depth}"
        raise ValueError(msg)
    # phyloseq drops samples with sample_sums(physeq) < sample.size; a sample at exactly depth stays.
    kept = np.flatnonzero(sums >= depth)
    if kept.size == 0:
        msg = f"no sample has depth={depth} reads; the deepest has {sums.max():g}"
        raise ValueError(msg)
    if kept.size < adata.n_obs:
        dropped = adata.obs_names[sums < depth].tolist()
        warn_user(f"pp.rarefy dropped {len(dropped)} sample(s) with fewer than depth={depth} reads: {dropped[:5]}")
    counts = _subsample_rows(X[kept], depth, as_generator(seed))
    present = np.unique(counts.indices).astype(np.intp)
    out = feature_subset(adata[kept], present)
    out.X = counts[:, present]
    add_provenance(out, "pp.rarefy", depth=depth)
    return out


def _smallest_nonzero(sums: npt.NDArray[np.int64]) -> int:
    # phyloseq's default is min(sample_sums(physeq)); skipping empty samples drops them instead of
    # rarefying everything to 0. With every sample empty this is 1, and the caller's no-sample check raises.
    nonzero = sums[sums > 0]
    return int(nonzero.min()) if nonzero.size else 1


def _subsample_rows(X: sp.csr_matrix, depth: int, rng: np.random.Generator) -> sp.csr_matrix:
    # Only the stored values are drawn from; zeros cannot be drawn, so positions carry over.
    data = np.concatenate(
        [
            subsample_counts(X.data[start:end].astype(np.int64), depth, replace=False, seed=rng)
            for start, end in zip(X.indptr[:-1], X.indptr[1:], strict=True)
        ]
    )
    # subsample_counts returns the platform's default int (int32 on Windows); pin int64 explicitly.
    out = sp.csr_matrix((data.astype(np.int64, copy=False), X.indices.copy(), X.indptr.copy()), shape=X.shape)
    out.eliminate_zeros()
    return out
