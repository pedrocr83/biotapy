from importlib.metadata import version

from . import datasets, pp

__all__ = ["__version__", "datasets", "pp"]

__version__ = version("biotapy")
