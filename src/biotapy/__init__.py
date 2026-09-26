from importlib.metadata import version

from . import datasets

__all__ = ["__version__", "datasets"]

__version__ = version("biotapy")
