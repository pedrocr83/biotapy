"""Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

from ._matrix import argmax_by, as_csr, sum_by
from ._optional import import_optional
from ._rng import as_generator
from ._slots import XKind, add_provenance, feature_subset, infer_x_kind, require_counts, x_kind
from ._taxonomy import RANKS, normalize_ranks, split_lineage, split_ranks
from ._tree import (
    PHYLO_KEY,
    TreeData,
    get_tree,
    make_treedata,
    relabel_tips,
    tree_from_edges,
    tree_from_newick,
    tree_tips,
)
from ._warnings import warn_user

__all__ = [
    "PHYLO_KEY",
    "RANKS",
    "TreeData",
    "XKind",
    "add_provenance",
    "argmax_by",
    "as_csr",
    "as_generator",
    "feature_subset",
    "get_tree",
    "import_optional",
    "infer_x_kind",
    "make_treedata",
    "normalize_ranks",
    "relabel_tips",
    "require_counts",
    "split_lineage",
    "split_ranks",
    "sum_by",
    "tree_from_edges",
    "tree_from_newick",
    "tree_tips",
    "warn_user",
    "x_kind",
]
