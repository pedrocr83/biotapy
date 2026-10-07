"""Differential abundance at 2,000 samples x 10,000 features: the two methods that run without R."""

from anndata import AnnData

import biotapy as bt

from ._data import synthetic_genera


class DifferentialAbundance:
    """``da.linda`` and ``da.ancombc2``, two groups of 1,000 samples; each densifies ``X`` once (160 MB)."""

    number, repeat, rounds, timeout = 1, 3, 1, 600

    def setup_cache(self) -> AnnData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic_genera()

    def time_linda(self, adata: AnnData) -> None:
        """One least-squares fit per feature, then the mode of the coefficients."""
        bt.da.linda(adata, "group")

    def peakmem_linda(self, adata: AnnData) -> None:
        """Peak memory; the docstring says about five dense copies of ``X``."""
        bt.da.linda(adata, "group")

    def time_ancombc2(self, adata: AnnData) -> None:
        """scikit-bio's ANCOM-BC2, its bias E-M capped at 100 iterations."""
        bt.da.ancombc2(adata, "group")

    def peakmem_ancombc2(self, adata: AnnData) -> None:
        """Peak memory; ``da.ancombc2``'s docstring gives 4.5x to 7.5x the dense ``X``, 4.5x on this table."""
        bt.da.ancombc2(adata, "group")
