import anndata as ad
import numpy as np
import pandas as pd
import pytest

from biotapy._core import split_ranks


def _adata(columns: list[str]) -> ad.AnnData:
    var = pd.DataFrame({c: ["x"] for c in columns}, index=["f1"])
    return ad.AnnData(X=np.zeros((1, 1)), obs=pd.DataFrame(index=["s1"]), var=var)


def test_split_ranks_uses_canonical_order():
    assert split_ranks(_adata(["genus", "kingdom", "phylum"]), "phylum") == (["kingdom", "phylum"], ["genus"])


def test_split_ranks_ignores_non_rank_columns():
    assert split_ranks(_adata(["kingdom", "sequence"]), "kingdom") == (["kingdom"], [])


def test_split_ranks_names_the_missing_rank():
    with pytest.raises(KeyError, match="genus"):
        split_ranks(_adata(["kingdom"]), "genus")
