"""x_kind, provenance and feature-changing subsets (contracts/data-model-slots)."""

import json
from importlib.metadata import version
from typing import Literal, cast

import numpy as np
import numpy.typing as npt
from anndata import AnnData

from ._matrix import as_csr

XKind = Literal["counts", "relative", "rpk", "cpm", "abundance"]
ParamValue = str | int | float | bool | None
DERIVED_SLOTS = ("layers", "obsm", "obsp", "varm", "varp")
# The uns['biotapy'] keys that survive a feature change; the others (pcoa, nmds) described the old features.
KEPT_META = ("x_kind", "provenance")
# Stored proportions are rounded: enterotype's sample sums run 0.99986-1.00000
# (decisions/phyloseq-import-route), so an exact test for 1 would call them abundances.
RELATIVE_TOLERANCE = 1e-3


def x_kind(adata: AnnData) -> XKind:
    """What ``X`` holds; a missing key means raw counts."""
    return cast(XKind, adata.uns.get("biotapy", {}).get("x_kind", "counts"))


def infer_x_kind(X: object) -> XKind:
    """What a freshly read ``X`` holds, judged from its values.

    Whole numbers are ``"counts"``; otherwise, rows that each sum to 1 within
    ``RELATIVE_TOLERANCE`` (all-zero rows ignored) are ``"relative"``;
    anything else is ``"abundance"``.
    """
    matrix = as_csr(X)
    if np.all(matrix.data == np.round(matrix.data)):
        return "counts"
    sums = np.asarray(matrix.sum(axis=1)).ravel()
    if np.all(np.abs(sums[sums != 0] - 1) <= RELATIVE_TOLERANCE):
        return "relative"
    return "abundance"


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
        # anndata 0.13 lists X itself as layers[None]; deleting that key would delete X.
        for key in [key for key in mapping.keys() if key is not None]:
            del mapping[key]
    meta = out.uns.get("biotapy", {"x_kind": "counts"})
    out.uns = {"biotapy": {key: meta[key] for key in KEPT_META if key in meta}}
    return out
