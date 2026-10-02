from importlib.metadata import version

from . import datasets, io, pl, pp, tl

__all__ = ["__version__", "datasets", "io", "pl", "pp", "tl"]

__version__ = version("biotapy")
