"""Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

from ._composition import check_pseudocount, pseudocounted
from ._download import make_pooch
from ._function import (
    BY_TAXON_KEY,
    FUNCTION_KEY,
    PROTECTED_FEATURES,
    SPECIAL_FEATURES,
    UNGROUPED,
    function_var,
    make_function_mudata,
)
from ._matrix import argmax_by, as_csr, divide_rows, finite_non_negative, sum_by, sum_pairs
from ._optional import import_optional
from ._rng import as_generator
from ._slots import (
    RELATIVE_TOLERANCE,
    XKind,
    add_provenance,
    feature_subset,
    infer_x_kind,
    replace_features,
    require_categorical,
    require_counts,
    x_kind,
)
from ._taxonomy import RANKS, normalize_ranks, split_lineage, split_ranks
from ._tree import (
    PHYLO_KEY,
    TreeData,
    get_skbio_tree,
    get_tree,
    make_treedata,
    relabel_tips,
    tree_from_edges,
    tree_from_newick,
    tree_from_phylo,
    tree_tips,
)
from ._warnings import warn_user

__all__ = [
    "BY_TAXON_KEY",
    "FUNCTION_KEY",
    "PHYLO_KEY",
    "PROTECTED_FEATURES",
    "RANKS",
    "RELATIVE_TOLERANCE",
    "SPECIAL_FEATURES",
    "UNGROUPED",
    "TreeData",
    "XKind",
    "add_provenance",
    "argmax_by",
    "as_csr",
    "as_generator",
    "check_pseudocount",
    "divide_rows",
    "feature_subset",
    "finite_non_negative",
    "function_var",
    "get_skbio_tree",
    "get_tree",
    "import_optional",
    "infer_x_kind",
    "make_function_mudata",
    "make_pooch",
    "make_treedata",
    "normalize_ranks",
    "pseudocounted",
    "relabel_tips",
    "replace_features",
    "require_categorical",
    "require_counts",
    "split_lineage",
    "split_ranks",
    "sum_by",
    "sum_pairs",
    "tree_from_edges",
    "tree_from_newick",
    "tree_from_phylo",
    "tree_tips",
    "warn_user",
    "x_kind",
]
