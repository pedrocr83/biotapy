"""Lazy imports for optional extras (decisions/optional-heavy-dependencies)."""

import importlib
from types import ModuleType


def import_optional(name: str, *, extra: str) -> ModuleType:
    """Import ``name`` or raise an ImportError that names the extra to install."""
    try:
        return importlib.import_module(name)
    except ImportError as err:
        msg = f"{name} is required here. Install it with: pip install 'biotapy[{extra}]'"
        raise ImportError(msg) from err
