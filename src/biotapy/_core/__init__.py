"""Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

from ._matrix import argmax_by, as_csr, sum_by
from ._optional import import_optional
from ._rng import as_generator
from ._slots import XKind, add_provenance, feature_subset, require_counts, x_kind
from ._taxonomy import RANKS, split_ranks

__all__ = [
    "RANKS",
    "XKind",
    "add_provenance",
    "argmax_by",
    "as_csr",
    "as_generator",
    "feature_subset",
    "import_optional",
    "require_counts",
    "split_ranks",
    "sum_by",
    "x_kind",
]
