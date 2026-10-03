"""Diversity at 5,000 samples x 50,000 features, and Bray-Curtis on sample subsets."""

import biotapy as bt
from biotapy._core import TreeData

from ._data import synthetic


class Alpha:
    """``tl.alpha``: rows densified in chunks of 2**20 values."""

    number, repeat, rounds, timeout = 1, 1, 1, 900

    def setup_cache(self) -> TreeData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic()

    def time_alpha(self, tdata: TreeData) -> None:
        """The four default metrics."""
        bt.tl.alpha(tdata)

    def peakmem_alpha(self, tdata: TreeData) -> None:
        """Peak memory of the four default metrics: chunking keeps it far below one dense copy (2 GB)."""
        bt.tl.alpha(tdata)

    def time_faith_pd(self, tdata: TreeData) -> None:
        """Faith PD on a 100,000-node tree, which scikit-bio re-indexes for each chunk of 20 samples."""
        bt.tl.alpha(tdata, metrics=["faith_pd"])


class Beta:
    """``tl.beta`` Bray-Curtis: pairwise, so time grows with the square of the samples."""

    number, repeat, rounds, timeout = 1, 1, 1, 1800
    params, param_names = [500, 1_000, 2_000], ["n_obs"]

    def setup_cache(self) -> TreeData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic()

    def time_beta(self, tdata: TreeData, n_obs: int) -> None:
        """The first ``n_obs`` samples, all 50,000 features."""
        bt.tl.beta(tdata[:n_obs].copy())


class BetaFullSize:
    """Peak memory of ``tl.beta`` on all 5,000 samples: X densified once (2 GB) plus the distances."""

    number, repeat, rounds, timeout = 1, 1, 1, 1800

    def setup_cache(self) -> TreeData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic()

    def peakmem_beta(self, tdata: TreeData) -> None:
        """One call, about 6 minutes."""
        bt.tl.beta(tdata)
