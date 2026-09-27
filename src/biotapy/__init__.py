from importlib.metadata import version

from . import datasets, io, pp, tl

__all__ = ["__version__", "datasets", "io", "pp", "tl"]

__version__ = version("biotapy")
