"""Preprocessing at 5,000 samples x 50,000 features, and PhILR at 1,000 samples."""

import biotapy as bt
from biotapy._core import TreeData

from ._data import N_PHILR_OBS, synthetic


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


class Philr:
    """``pp.philr`` at 1,000 samples x 50,000 features: ``X`` densified once (400 MB), then 49,999 balances."""

    number, repeat, rounds, timeout = 1, 3, 1, 600

    def setup_cache(self) -> TreeData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic(n_obs=N_PHILR_OBS)

    def time_philr(self, tdata: TreeData) -> None:
        """Centred log-ratios times the tree's sparse orthonormal basis."""
        bt.pp.philr(tdata)

    def peakmem_philr(self, tdata: TreeData) -> None:
        """Peak memory: the dense table, its log-ratios and the balances."""
        bt.pp.philr(tdata)
