"""Renormalisation of HUMAnN tables, with humann_renorm_table's community semantics."""

from typing import Literal, cast

import numpy as np
from anndata import AnnData
from mudata import MuData

from biotapy._core import (
    BY_TAXON_KEY,
    FUNCTION_KEY,
    SPECIAL_FEATURES,
    add_provenance,
    as_csr,
    divide_rows,
    feature_subset,
    warn_user,
)

_SCALE = {"relab": 1.0, "cpm": 1e6}
_X_KIND: dict[str, Literal["relative", "cpm"]] = {"relab": "relative", "cpm": "cpm"}


def renorm(mdata: MuData, units: Literal["relab", "cpm"], *, special: bool = True) -> MuData:
    """Rescale both modalities so each sample's community total is 1 (or one million).

    Parameters
    ----------
    mdata
        ``bt.io.read_humann``'s result: modalities ``"function"`` and
        ``"function_by_taxon"``.
    units
        ``"relab"`` (totals 1) or ``"cpm"`` (totals 1,000,000).
    special
        Keep ``UNMAPPED``, ``READS_UNMAPPED``, ``UNINTEGRATED`` and
        ``UNGROUPED`` (community and per-taxon rows) and count them in the
        total; ``False`` drops them first.

    Returns
    -------
    MuData
        A copy whose two function modalities have ``X`` divided by each
        sample's total over the ``"function"`` modality, and
        ``uns['biotapy']['x_kind']`` set to ``"relative"`` or ``"cpm"``.
        Their ``layers``, ``obsm``, ``obsp``, ``varm`` and ``varp`` are
        dropped, as after any change to ``X``'s features. Other modalities
        and the global ``obs`` are copied unchanged.

    Raises
    ------
    KeyError
        ``mdata`` lacks one of the two modalities.
    ValueError
        ``units`` is not ``"relab"`` or ``"cpm"``; the two modalities do not
        hold the same samples in the same order; the ``"function"``
        modality has no feature, or none left after ``special=False``.

    Warns
    -----
    UserWarning
        Some samples have a zero community total; the warning names up to
        three of them and the count. They stay zero, as in HUMAnN.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/function`

    Matches ``humann_renorm_table`` (HUMAnN 3.9) in its default community
    mode: stratified rows are divided by the community total, so a
    pathway's strata need not sum to its community value. For HUMAnN's
    ``--mode levelwise``, use ``bt.pp.relative`` on each modality. A sample
    with a zero total stays zero.

    References
    ----------
    Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of
    diverse microbial communities with bioBakery 3. eLife 10:e65088.

    Examples
    --------
    >>> import biotapy as bt
    >>> out = bt.fn.renorm(bt.datasets.toy_humann(), "relab")
    >>> round(float(out["function"].X[0].sum()), 6)
    1.0
    """
    if units not in _SCALE:
        msg = f"units={units!r} must be 'relab' or 'cpm'"
        raise ValueError(msg)
    community, by_taxon = (_kept(mod, special=special) for mod in _function_modalities(mdata))
    if community.n_vars == 0:
        reason = (
            "special=False dropped every community row, as all of them are special"
            if not special and mdata.mod[FUNCTION_KEY].n_vars > 0
            else "renorm needs the community (unstratified) rows"
        )
        msg = f"mdata[{FUNCTION_KEY!r}] has no feature to total; {reason}"
        raise ValueError(msg)
    totals = np.asarray(as_csr(community.X).sum(axis=1), dtype=np.float64).ravel()
    empty = community.obs_names[totals == 0]
    if len(empty) > 0:
        warn_user(f"fn.renorm: {len(empty)} sample(s) have no community abundance and stay zero: {list(empty[:3])}")
    out = mdata.copy()
    # ModDict is a dict; mudata types the property as a read-only Mapping.
    mods = cast("dict[str, AnnData | MuData]", out.mod)
    for key, mod in ((FUNCTION_KEY, community), (BY_TAXON_KEY, by_taxon)):
        mods[key] = _rescaled(mod, totals, units=units, special=special)
    out.update()
    return out


def _function_modalities(mdata: MuData) -> tuple[AnnData, AnnData]:
    community, by_taxon = mdata.mod.get(FUNCTION_KEY), mdata.mod.get(BY_TAXON_KEY)
    if not isinstance(community, AnnData) or not isinstance(by_taxon, AnnData):
        msg = f"mdata needs AnnData modalities {[FUNCTION_KEY, BY_TAXON_KEY]}, as bt.io.read_humann returns them"
        raise KeyError(msg)
    if not community.obs_names.equals(by_taxon.obs_names):
        msg = (
            f"mdata[{FUNCTION_KEY!r}] and mdata[{BY_TAXON_KEY!r}] must hold the same samples in the same order, "
            "because each sample's community total divides its own stratified rows"
        )
        raise ValueError(msg)
    return community, by_taxon


def _kept(adata: AnnData, *, special: bool) -> AnnData:
    """``adata`` without its special rows unless ``special``; a real object, never a view."""
    function = adata.var["function"] if "function" in adata.var.columns else adata.var_names.to_series()
    keep = np.ones(adata.n_vars, dtype=bool) if special else ~function.isin(SPECIAL_FEATURES).to_numpy()
    return feature_subset(adata, np.flatnonzero(keep))


def _rescaled(adata: AnnData, totals: np.ndarray, *, units: str, special: bool) -> AnnData:
    """Divide each stored value by its sample's total, then scale to ``units``."""
    adata.X = divide_rows(as_csr(adata.X), totals) * _SCALE[units]
    adata.uns["biotapy"]["x_kind"] = _X_KIND[units]
    add_provenance(adata, "fn.renorm", units=units, special=special)
    return adata
