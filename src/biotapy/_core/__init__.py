"""Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

from ._optional import import_optional
from ._rng import as_generator

__all__ = ["as_generator", "import_optional"]
