"""Preprocessing at 5,000 samples x 50,000 features."""

import biotapy as bt
from biotapy._core import TreeData

from ._data import synthetic


class Preprocessing:
    """``relative``, ``tax_glom`` and ``rarefy`` on the full synthetic table."""

    number, repeat, rounds, timeout = 1, 3, 1, 300

    def setup_cache(self) -> TreeData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic()

    def time_relative(self, tdata: TreeData) -> None:
        """``pp.relative``."""
        bt.pp.relative(tdata)

    def time_tax_glom(self, tdata: TreeData) -> None:
        """``pp.tax_glom`` to 1,000 genera."""
        bt.pp.tax_glom(tdata, "genus")

    def time_rarefy(self, tdata: TreeData) -> None:
        """``pp.rarefy`` to the smallest sample depth."""
        bt.pp.rarefy(tdata, seed=0)
