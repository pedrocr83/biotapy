"""x_kind, provenance and feature-changing subsets (contracts/data-model-slots)."""

import json
from importlib.metadata import version
from typing import Literal, cast

import numpy as np
import numpy.typing as npt
from anndata import AnnData

XKind = Literal["counts", "relative", "rpk", "cpm", "abundance"]
ParamValue = str | int | float | bool | None
DERIVED_SLOTS = ("layers", "obsm", "obsp", "varm", "varp")


def x_kind(adata: AnnData) -> XKind:
    """What ``X`` holds; a missing key means raw counts."""
    return cast(XKind, adata.uns.get("biotapy", {}).get("x_kind", "counts"))


def require_counts(adata: AnnData, *, func: str) -> None:
    """Raise unless ``X`` holds raw counts."""
    kind = x_kind(adata)
    if kind != "counts":
        msg = f"{func} needs raw counts in X, but uns['biotapy']['x_kind'] is {kind!r}"
        raise ValueError(msg)


def add_provenance(adata: AnnData, step: str, **params: ParamValue) -> None:
    """Append one JSON entry to ``uns['biotapy']['provenance']`` (h5ad cannot store a list of dicts)."""
    meta = adata.uns.setdefault("biotapy", {"x_kind": "counts"})
    entry = json.dumps({"step": step, "version": version("biotapy"), "params": params})
    meta["provenance"] = [*meta.get("provenance", []), entry]


def feature_subset(adata: AnnData, index: npt.NDArray[np.intp]) -> AnnData:
    """Subset features and drop every slot derived from the old feature set."""
    out = adata[:, index].copy()
    for slot in DERIVED_SLOTS:
        mapping = getattr(out, slot)
        for key in list(mapping.keys()):
            del mapping[key]
    out.uns = {"biotapy": out.uns.get("biotapy", {"x_kind": "counts"})}
    return out
