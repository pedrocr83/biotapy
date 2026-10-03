from importlib.metadata import version

from . import datasets, fn, io, pl, pp, tl

__all__ = ["__version__", "datasets", "fn", "io", "pl", "pp", "tl"]

__version__ = version("biotapy")
