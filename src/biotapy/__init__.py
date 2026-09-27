from importlib.metadata import version

from . import datasets, io, pp

__all__ = ["__version__", "datasets", "io", "pp"]

__version__ = version("biotapy")
