from importlib.metadata import version

from . import da, datasets, fn, io, pl, pp, tl

__all__ = ["__version__", "da", "datasets", "fn", "io", "pl", "pp", "tl"]

__version__ = version("biotapy")
