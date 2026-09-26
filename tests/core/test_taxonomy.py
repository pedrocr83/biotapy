import anndata as ad
import numpy as np
import pandas as pd
import pytest

from biotapy._core import normalize_ranks, split_lineage, split_ranks


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


LINEAGES = pd.Series(
    [
        "k__Bacteria; p__Firmicutes; c__; o__; f__; g__; s__",
        "d__Bacteria; p__Proteobacteria",
        "D_0__Archaea;D_1__Euryarchaeota",
        "Bacteria;Firmicutes",
        "p__Firmicutes_A;c__Clostridia_258483",
        np.nan,
        "Unassigned",
    ],
    index=[f"f{i}" for i in range(7)],
)


def test_split_lineage_reads_greengenes_prefixes():
    out = split_lineage(LINEAGES)
    assert list(out.columns) == ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    assert out.loc["f0", ["kingdom", "phylum"]].tolist() == ["Bacteria", "Firmicutes"]
    assert out.loc["f0", "class":].isna().all()


@pytest.mark.parametrize(
    ("feature", "expected"),
    [("f1", ["Bacteria", "Proteobacteria"]), ("f2", ["Archaea", "Euryarchaeota"]), ("f3", ["Bacteria", "Firmicutes"])],
)
def test_split_lineage_reads_silva_and_unprefixed(feature, expected):
    assert split_lineage(LINEAGES).loc[feature, ["kingdom", "phylum"]].tolist() == expected


def test_split_lineage_places_truncated_lineage_by_prefix():
    row = split_lineage(LINEAGES).loc["f4"]
    assert pd.isna(row["kingdom"]) and row["phylum"] == "Firmicutes_A" and row["class"] == "Clostridia_258483"


def test_split_lineage_missing_lineage_is_all_nan():
    assert split_lineage(LINEAGES).loc["f5"].isna().all()


def test_split_lineage_keeps_unassigned_as_kingdom():
    assert split_lineage(LINEAGES).loc["f6", "kingdom"] == "Unassigned"


def test_split_lineage_columns_stop_at_deepest_rank():
    assert list(split_lineage(pd.Series(["k__A; p__B"], index=["x"])).columns) == ["kingdom", "phylum"]


def test_normalize_ranks_canonicalizes_names_and_missing_values():
    frame = pd.DataFrame(
        {"Domain": ["Bacteria", "NA", " "], "Genus": ["g__", None, "g__Blautia"], "sequence": ["AC", "GT", "TT"]},
        index=["a", "b", "c"],
    )
    out = normalize_ranks(frame)
    assert list(out.columns) == ["kingdom", "genus", "sequence"]
    assert out.loc["a", "kingdom"] == "Bacteria" and out["kingdom"].iloc[1:].isna().all()
    assert out["genus"].iloc[:2].isna().all() and out.loc["c", "genus"] == "Blautia"
    assert out["sequence"].tolist() == ["AC", "GT", "TT"]


def test_normalize_ranks_leaves_input_alone():
    frame = pd.DataFrame({"Genus": ["g__Blautia"]})
    normalize_ranks(frame)
    assert list(frame.columns) == ["Genus"] and frame.iloc[0, 0] == "g__Blautia"
