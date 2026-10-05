"""Function tables: reading HUMAnN output, regrouping it, and Tian et al.'s functional redundancy."""

from pathlib import Path

import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

import biotapy as bt

from ._data import function_groups, synthetic_function, synthetic_traits


class FuncGlom:
    """``fn.func_glom`` of a 1,600 x 22,000 pathway table, every function in two of 50 groups."""

    number, repeat, rounds, timeout = 1, 3, 1, 300

    def setup_cache(self) -> dict[str, AnnData]:
        """Build the two modalities once; asv pickles them for every benchmark (a MuData does not unpickle)."""
        return dict(synthetic_function().mod)

    def setup(self, modalities: dict[str, AnnData]) -> None:
        """The edge table is small; build it outside the timed call."""
        self.groups = function_groups()

    def time_func_glom_community(self, modalities: dict[str, AnnData]) -> None:
        """The 500 community rows to 50 groups."""
        bt.fn.func_glom(modalities["function"], "group", hierarchy=self.groups)

    def time_func_glom_by_taxon(self, modalities: dict[str, AnnData]) -> None:
        """The 21,500 stratified rows to 50 groups per taxon."""
        bt.fn.func_glom(modalities["function_by_taxon"], "group", hierarchy=self.groups)


class ReadHumann:
    """``io.read_humann`` of the same table written as HUMAnN writes a merged pathway table."""

    number, repeat, rounds, timeout = 1, 3, 1, 600

    def setup_cache(self) -> str:
        """Write the table once, in asv's working directory for this class."""
        mdata = synthetic_function()
        values = sp.hstack([mdata["function"].X, mdata["function_by_taxon"].X], format="csc").T.toarray()
        ids = [*mdata["function"].var_names, *mdata["function_by_taxon"].var_names]
        table = pd.DataFrame(
            values, index=pd.Index(ids, name="# Pathway"), columns=[f"{s}_Abundance" for s in mdata.obs_names]
        )
        path = Path("pathabundance.tsv").resolve()
        table.to_csv(path, sep="\t", float_format="%.6g")
        return str(path)

    def time_read_humann(self, path: str) -> None:
        """22,000 rows x 1,600 samples."""
        bt.io.read_humann(path)

    def peakmem_read_humann(self, path: str) -> None:
        """Peak memory: the reader builds one dense rows x samples float64 array (282 MB) before CSR."""
        bt.io.read_humann(path)


class FunctionalRedundancy:
    """``fn.functional_redundancy`` at 2,000 taxa: pairwise taxon distances, so O(taxa^2) time and memory."""

    number, repeat, rounds, timeout = 1, 1, 1, 600

    def setup_cache(self) -> tuple[AnnData, pd.DataFrame]:
        """Build the abundances and the copy numbers once."""
        return synthetic_traits()

    def time_functional_redundancy(self, data: tuple[AnnData, pd.DataFrame]) -> None:
        """100 samples, 2,000 taxa, 2,500 genes."""
        adata, traits = data
        bt.fn.functional_redundancy(adata, traits=traits)

    def peakmem_functional_redundancy(self, data: tuple[AnnData, pd.DataFrame]) -> None:
        """Peak memory, dominated by the 2,000 x 2,000 distance matrix."""
        adata, traits = data
        bt.fn.functional_redundancy(adata, traits=traits)
