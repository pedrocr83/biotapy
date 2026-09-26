"""The only module that imports treedata or networkx (contracts/tree-access)."""

# networkx's real DiGraph class is not subscriptable at runtime (types-networkx is a
# stub-only package); postponed evaluation keeps the `nx.DiGraph[str]` annotations
# below from being evaluated eagerly and raising TypeError.
from __future__ import annotations

from collections.abc import Iterable
from typing import cast

import networkx as nx
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from treedata import TreeData as TreeData

from ._matrix import as_csr
from ._slots import XKind, add_provenance

PHYLO_KEY = "phylo"


def tree_from_edges(edges: Iterable[tuple[str, str, float]]) -> nx.DiGraph[str]:
    """Build a rooted tree from ``(parent, child, branch_length)`` triples."""
    tree: nx.DiGraph[str] = nx.DiGraph()
    tree.add_weighted_edges_from(edges, weight="length")
    return tree


def get_tree(tdata: TreeData) -> nx.DiGraph[str]:
    """Return the phylogeny in ``vart['phylo']``."""
    if PHYLO_KEY not in tdata.vart:
        msg = f"no phylogeny in vart[{PHYLO_KEY!r}]"
        raise KeyError(msg)
    # treedata ships no py.typed marker, so tdata.vart[...] types as Any.
    return cast("nx.DiGraph[str]", tdata.vart[PHYLO_KEY])


def make_treedata(
    X: sp.spmatrix | npt.ArrayLike,
    *,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    tree: nx.DiGraph[str] | None,
    x_kind: XKind,
    source: str,
) -> TreeData:
    """Construct a TreeData that follows contracts/data-model-slots."""
    vart = None if tree is None else {PHYLO_KEY: tree}
    tdata = TreeData(X=as_csr(X), obs=obs, var=var, vart=vart, label=None)
    tdata.uns["biotapy"] = {"x_kind": x_kind}
    add_provenance(tdata, source)
    return tdata
