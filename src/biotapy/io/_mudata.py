"""Several data types over the same samples as one MuData (decisions/multiomics-as-mudata)."""

from collections.abc import Mapping

import pandas as pd
from anndata import AnnData
from mudata import MuData

from biotapy._core import warn_user


def to_mudata(modalities: Mapping[str, AnnData]) -> MuData:
    """Combine data types measured on the same samples into one MuData.

    Parameters
    ----------
    modalities
        Modality name -> samples x features table. biotapy's names are
        ``"taxa"``, ``"function"``, ``"function_by_taxon"``, ``"metabolites"``
        and ``"host"``; a TreeData stays a TreeData.

    Returns
    -------
    MuData
        One modality per entry, in the mapping's order, each a copy holding only
        the samples every modality has, in the first modality's order. The
        global ``obs`` has no columns; ``mdata.pull_obs()`` gathers them.

    Raises
    ------
    TypeError
        A name is not a non-empty string, or a value is not an AnnData (such as
        a function table, which is already a MuData).
    ValueError
        ``modalities`` is empty, a modality repeats a sample id, or the
        modalities share no sample.

    Warns
    -----
    UserWarning
        Some samples are missing from at least one modality; the warning names
        how many each modality loses.

    Notes
    -----
    R equivalent: ``MultiAssayExperiment::MultiAssayExperiment``
    Guide: :doc:`/guide/multiomics`

    A function table from ``bt.io.read_humann`` or ``bt.io.read_picrust2`` is
    a MuData of two modalities; pass them as two entries, ``{**table.mod,
    "taxa": tdata}``.

    ``MuData.write_h5mu`` does not keep a TreeData's tree: the modality reads
    back as an AnnData without ``vart``. Save a tree-bearing modality with
    ``TreeData.write_h5td`` as well.

    Examples
    --------
    >>> import biotapy as bt
    >>> table = bt.datasets.toy_humann()
    >>> mdata = bt.io.to_mudata({**table.mod, "taxa": bt.datasets.toy()})
    >>> list(mdata.mod)
    ['function', 'function_by_taxon', 'taxa']
    """
    if not modalities:
        msg = "modalities is empty; pass at least one table, e.g. {'taxa': tdata}"
        raise ValueError(msg)
    shared: pd.Index[str] | None = None
    for name, mod in modalities.items():
        _check_modality(name, mod)
        shared = mod.obs_names if shared is None else shared.intersection(mod.obs_names, sort=False)
    if shared is None or shared.empty:
        firsts = {name: mod.obs_names[:3].tolist() for name, mod in modalities.items()}
        msg = f"modalities share no sample; their first sample ids: {firsts}"
        raise ValueError(msg)
    lost = [
        f"{name} {mod.n_obs - len(shared)} of {mod.n_obs}"
        for name, mod in modalities.items()
        if mod.n_obs > len(shared)
    ]
    if lost:
        warn_user(f"samples missing from another modality are dropped: {', '.join(lost)}")
    return MuData({name: mod[shared].copy() for name, mod in modalities.items()})


def _check_modality(name: object, mod: object) -> None:
    """Raise unless ``name`` is a non-empty string and ``mod`` an AnnData with unique sample ids."""
    if not isinstance(name, str) or not name:
        msg = f"modality names must be non-empty strings, got {name!r}"
        raise TypeError(msg)
    if not isinstance(mod, AnnData):
        msg = (
            f"modalities[{name!r}] must be an AnnData, got {type(mod).__name__}; "
            "pass a function table's modalities as two entries: {**table.mod, ...}"
        )
        raise TypeError(msg)
    repeated = mod.obs_names[mod.obs_names.duplicated()].unique().tolist()
    if repeated:
        msg = f"modalities[{name!r}] repeats sample ids: {repeated[:5]}"
        raise ValueError(msg)
