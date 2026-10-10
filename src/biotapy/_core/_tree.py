"""The only module that imports treedata or networkx (contracts/tree-access)."""

# networkx's real DiGraph class is not subscriptable at runtime (types-networkx is a
# stub-only package); postponed evaluation keeps the `nx.DiGraph[str]` annotations
# below from being evaluated eagerly and raising TypeError.
from __future__ import annotations

import io
import itertools
import math
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from typing import cast

import h5py
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from anndata.io import write_elem
from skbio import TreeNode
from skbio.io import NewickFormatError, UnrecognizedFormatError
from treedata import TreeData as TreeData
from treedata import read_h5td

from ._matrix import as_csr
from ._slots import XKind, add_provenance
from ._warnings import warn_user

PHYLO_KEY = "phylo"
TREE_SLOTS_ATTR = "biotapy-treedata-encoding"
_TREE_SLOTS_VERSION = "1"
_TREE_SLOTS = ("obst", "vart", "label", "allow_overlap", "alignment")
# treedata's root attributes, kept under biotapy names so the copied elements are read back as written.
_ROOT_ATTRS = ("encoding-type", "encoding-version")


def tree_from_edges(edges: Iterable[tuple[str, str, float]]) -> nx.DiGraph[str]:
    """Build a rooted tree from ``(parent, child, branch_length)`` triples."""
    tree: nx.DiGraph[str] = nx.DiGraph()
    tree.add_weighted_edges_from(edges, weight="length")
    return tree


def tree_tips(tree: nx.DiGraph[str]) -> list[str]:
    """Names of the tree's tips (nodes without children)."""
    return [node for node in tree.nodes if tree.out_degree(node) == 0]


def relabel_tips(tree: nx.DiGraph[str], names: Mapping[str, str]) -> nx.DiGraph[str]:
    """Rename the nodes listed in ``names`` (e.g. sequence -> ASV id); others keep theirs.

    A new name that already names a node which is not itself renamed raises
    ``ValueError``: networkx would silently merge the two nodes.
    """
    renamed = {old: new for old, new in names.items() if old in tree}
    clashes = [new for new in renamed.values() if new in tree and new not in renamed]
    if clashes:
        msg = f"names= would merge nodes: {clashes[:3]} already name other nodes of the tree"
        raise ValueError(msg)
    # types-networkx's overloads for relabel_nodes resolve to Any here (PEP 696 default type
    # params on a generic base class), even though the runtime call is exactly this typed.
    return cast("nx.DiGraph[str]", nx.relabel_nodes(tree, renamed, copy=True))


def get_tree(tdata: TreeData) -> nx.DiGraph[str]:
    """Return the phylogeny in ``vart['phylo']``."""
    if PHYLO_KEY not in tdata.vart:
        msg = f"no phylogeny in vart[{PHYLO_KEY!r}]"
        raise KeyError(msg)
    # treedata ships no py.typed marker, so tdata.vart[...] types as Any.
    return cast("nx.DiGraph[str]", tdata.vart[PHYLO_KEY])


def get_skbio_tree(adata: AnnData, *, split_root: bool = True) -> TreeNode:
    """The phylogeny in ``vart['phylo']`` as a scikit-bio ``TreeNode``, rooted where it is drawn.

    scikit-bio's Faith PD and UniFrac accept a root with at most two children. A
    root with more (an unrooted Newick tree, or the toy tree) keeps its first child
    and gets the others under one new zero-length node, which changes no
    root-to-tip path length; ``split_root=False`` leaves it as stored (PhILR must
    refuse it). Children keep the order of the stored edges. Missing (NaN) branch
    lengths stay NaN; scikit-bio counts them as zero. A plain AnnData raises
    ``TypeError``; a TreeData without ``vart['phylo']`` raises ``KeyError``.
    """
    if not isinstance(adata, TreeData):
        msg = f"needs a TreeData with a tree in vart[{PHYLO_KEY!r}], got {type(adata).__name__}"
        raise TypeError(msg)
    tree = get_tree(adata)
    root = next(node for node in tree if tree.in_degree(node) == 0)
    nodes = {root: TreeNode(name=root)}
    for parent, child in nx.bfs_edges(tree, root):
        nodes[child] = TreeNode(name=child, length=tree.edges[parent, child]["length"])
        nodes[parent].append(nodes[child])
    top = nodes[root]
    if split_root and len(top.children) > 2:
        split = TreeNode(length=0.0)
        split.extend(top.children[1:])
        top.append(split)
    return top


def tree_from_newick(text: str, *, argument: str = "text") -> nx.DiGraph[str]:
    """Parse one Newick tree; tips keep their names, internal nodes get unique ones.

    Internal labels (often support values such as ``0.95``) repeat, so they
    cannot name graph nodes and are dropped. A missing branch length is NaN.
    Malformed text, or a tip without a unique name, raises ``ValueError``
    naming ``argument`` (a reader passes e.g. ``"tree='tree.nwk'"``).
    """
    try:
        # skbio turns "_" into " " in unquoted names by default; ASV ids need them intact.
        root = TreeNode.read([text], convert_underscores=False)
    except (NewickFormatError, UnrecognizedFormatError) as error:
        # skbio's format sniffer rejects most malformed text (UnrecognizedFormatError)
        # before its Newick parser can raise NewickFormatError.
        msg = f"{argument} is not a valid Newick tree"
        raise ValueError(msg) from error
    tips = [tip.name for tip in root.tips()]
    _require_unique_names(tips, argument)
    names = _node_names(root, set(tips))
    tree = tree_from_edges(
        (names[id(node.parent)], names[id(node)], math.nan if node.length is None else float(node.length))
        for node in root.preorder(include_self=False)
    )
    # A lone root that is also a tip has no edges; add it explicitly so it still becomes a node.
    tree.add_nodes_from(names.values())
    return tree


def _require_unique_names(tips: list[str | None], argument: str) -> None:
    bad = [name for name, count in Counter(tips).items() if name is None or count > 1]
    if bad:
        msg = f"{argument} needs unique tip names; unnamed or repeated: {bad[:5]}"
        raise ValueError(msg)


def _node_names(root: TreeNode, tips: set[str]) -> dict[int, str]:
    fresh = (name for name in (f"n{i}" for i in itertools.count()) if name not in tips)
    # _require_unique_names already rejected None/duplicate tip names before this runs.
    return {id(node): cast(str, node.name) if node.is_tip() else next(fresh) for node in root.preorder()}


def tree_from_phylo(
    edge: npt.ArrayLike, lengths: npt.ArrayLike | None, tips: Sequence[str], *, argument: str = "tips"
) -> nx.DiGraph[str]:
    """Build a tree from an ape ``phylo`` edge matrix (1-based; tips are ``1..len(tips)``).

    Internal nodes get the same collision-free ``n<i>`` names as ``tree_from_newick``;
    missing branch lengths are NaN.
    """
    edges = np.asarray(edge, dtype=np.int64).reshape(-1, 2)
    names = [str(tip) for tip in tips]
    # list[str] is a Sequence[str | None] at runtime; invariant List typing just can't see it.
    _require_unique_names(cast("list[str | None]", names), argument)
    taken = set(names)
    fresh = (name for name in (f"n{i}" for i in itertools.count()) if name not in taken)
    internal = {int(node): next(fresh) for node in np.unique(edges) if node > len(names)}
    label = {**{i + 1: tip for i, tip in enumerate(names)}, **internal}
    length = np.full(len(edges), np.nan) if lengths is None else np.asarray(lengths, dtype=np.float64)
    tree = tree_from_edges(
        (label[int(parent)], label[int(child)], float(value))
        for (parent, child), value in zip(edges, length, strict=True)
    )
    tree.add_nodes_from(names)
    return tree


def make_treedata(
    X: object,
    *,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    tree: nx.DiGraph[str] | None,
    x_kind: XKind,
    source: str,
) -> TreeData:
    """Construct a TreeData that follows contracts/data-model-slots.

    Ids become unique strings. With a tree, only features that are its tips
    are kept and tips outside the table are pruned; a mismatch warns once.
    """
    obs, var = _with_str_ids(obs, "obs"), _with_str_ids(var, "var")
    matrix = as_csr(X)
    if tree is not None:
        matrix, var, tree = _align_tree(matrix, var, tree)
    vart = None if tree is None else {PHYLO_KEY: tree}
    tdata = TreeData(X=matrix, obs=obs, var=var, vart=vart, label=None)
    tdata.uns["biotapy"] = {"x_kind": x_kind}
    add_provenance(tdata, source)
    return tdata


def _with_str_ids(frame: pd.DataFrame, axis: str) -> pd.DataFrame:
    n_missing = int(np.count_nonzero(frame.index.isna()))
    if n_missing:
        msg = f"{axis} ids must not be missing; {n_missing} id(s) are NaN or None"
        raise ValueError(msg)
    out = frame.copy()
    out.index = out.index.astype(str)
    duplicated = out.index[out.index.duplicated()].unique().tolist()
    if duplicated:
        msg = f"duplicate {axis} ids: {duplicated[:5]}"
        raise ValueError(msg)
    return out


def _align_tree(
    X: sp.csr_matrix, var: pd.DataFrame, tree: nx.DiGraph[str]
) -> tuple[sp.csr_matrix, pd.DataFrame, nx.DiGraph[str]]:
    tips = set(tree_tips(tree))
    shared = var.index.isin(tips)
    n_extra = len(tips - set(var.index))
    if shared.all() and n_extra == 0:
        return X, var, tree
    if not shared.any():
        msg = (
            "no feature of the table is a tip of the tree; "
            f"features: {var.index[:3].tolist()}, tips: {sorted(tips)[:3]}"
        )
        raise ValueError(msg)
    msg = (
        f"tree and table disagree: {int((~shared).sum())} feature(s) not in the tree and "
        f"{n_extra} tree tip(s) not in the table; keeping the {int(shared.sum())} shared features"
    )
    warn_user(msg)
    kept = var.index[shared]
    keep_nodes = set(kept).union(*(nx.ancestors(tree, tip) for tip in kept))
    return X[:, np.flatnonzero(shared)], var.loc[kept], tree.subgraph(keep_nodes).copy()


def write_tree_slots(group: h5py.Group, tdata: TreeData) -> None:
    """Copy a TreeData's trees and tree settings into an open h5 group, marked with a layout version.

    treedata's own writer only writes a file root, so a slim copy (names, trees, ``label``,
    ``allow_overlap``, ``alignment``; no ``X``) goes through ``write_h5td`` into an in-memory
    file and its elements are copied over. ``tdata`` is not changed.
    """
    slim = TreeData(
        obs=pd.DataFrame(index=tdata.obs_names),
        var=pd.DataFrame(index=tdata.var_names),
        obst=dict(tdata.obst),
        vart=dict(tdata.vart),
        label=tdata.label,
        allow_overlap=tdata.allow_overlap,
        alignment=tdata.alignment,
    )
    buffer = io.BytesIO()
    # treedata annotates filename as a path; h5py.File takes the in-memory buffer it forwards it to.
    slim.write_h5td(buffer)  # type: ignore[arg-type]
    with h5py.File(buffer, "r") as source:
        for name in _TREE_SLOTS:
            # A newer mudata/treedata may already have written this element; ours is the layout read_tree_slots reads.
            if name in group:
                del group[name]
            source.copy(name, group)
        for attr in _ROOT_ATTRS:
            group.attrs[f"biotapy-treedata-root-{attr}"] = source.attrs[attr]
    group.attrs[TREE_SLOTS_ATTR] = _TREE_SLOTS_VERSION


def read_tree_slots(group: h5py.Group, adata: AnnData) -> TreeData:
    """Rebuild a TreeData from ``adata`` and the tree slots ``write_tree_slots`` put in ``group``."""
    if TREE_SLOTS_ATTR not in group.attrs:
        msg = f"group {group.name!r} has no {TREE_SLOTS_ATTR!r} attribute, so it holds no tree slots"
        raise KeyError(msg)
    version = str(group.attrs[TREE_SLOTS_ATTR])
    if version != _TREE_SLOTS_VERSION:
        msg = (
            f"group {group.name!r} has tree-slot layout version {version!r}, but this biotapy reads version "
            f"{_TREE_SLOTS_VERSION!r}: the file was written by a newer biotapy, so upgrade biotapy"
        )
        raise ValueError(msg)
    buffer = io.BytesIO()
    with h5py.File(buffer, "w") as target:
        for attr in _ROOT_ATTRS:
            target.attrs[attr] = group.attrs[f"biotapy-treedata-root-{attr}"]
        for name in _TREE_SLOTS:
            group.copy(name, target)
        write_elem(target, "obs", pd.DataFrame(index=adata.obs_names))
        write_elem(target, "var", pd.DataFrame(index=adata.var_names))
    buffer.seek(0)
    slots = read_h5td(buffer)  # type: ignore[arg-type]
    # AnnData's constructor accepts an AnnData as X and copies every slot; treedata's stub narrows it.
    return TreeData(
        adata,  # type: ignore[arg-type]
        obst=dict(slots.obst),
        vart=dict(slots.vart),
        label=slots.label,
        allow_overlap=slots.allow_overlap,
        alignment=slots.alignment,
    )
