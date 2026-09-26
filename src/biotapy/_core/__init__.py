"""Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

from ._matrix import argmax_by, as_csr, sum_by
from ._optional import import_optional
from ._rng import as_generator

__all__ = ["argmax_by", "as_csr", "as_generator", "import_optional", "sum_by"]
