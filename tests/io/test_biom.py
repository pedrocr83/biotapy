import tempfile
from pathlib import Path

import anndata as ad
import biom
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import treedata
from biom.exception import TableException
from biom.util import biom_open
from hypothesis import given, settings
from hypothesis import strategies as st

import biotapy as bt

RANK_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus"]
ALL_RANKS = [*RANK_COLUMNS, "species"]


def test_read_biom_hdf5_puts_samples_in_rows(biom_hdf5, biom_table):
    tdata = bt.io.read_biom(biom_hdf5)
    assert tdata.shape == (3, 4)
    assert list(tdata.obs_names) == ["S1", "S2", "S3"]
    assert list(tdata.var_names) == ["OTU_1", "OTU_2", "OTU_3", "OTU_4"]
    np.testing.assert_array_equal(tdata.X.toarray(), biom_table.matrix_data.T.toarray())


def test_read_biom_json_matches_hdf5(biom_json, biom_hdf5):
    from_json, from_hdf5 = bt.io.read_biom(biom_json), bt.io.read_biom(biom_hdf5)
    assert (from_json.X != from_hdf5.X).nnz == 0
    pd.testing.assert_frame_equal(from_json.var, from_hdf5.var)


def test_read_biom_taxonomy_becomes_rank_columns(biom_hdf5):
    var = bt.io.read_biom(biom_hdf5).var
    assert list(var.columns) == ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    assert var.loc["OTU_2", "order"] == "Lactobacillales"
    assert var.loc["OTU_1", ["order", "family", "genus", "species"]].isna().all()


def test_read_biom_sample_metadata_goes_to_obs(biom_hdf5):
    assert bt.io.read_biom(biom_hdf5).obs["group"].tolist() == ["A", "A", "B"]


def test_read_biom_integer_ids_become_strings(biom_json_ids):
    tdata = bt.io.read_biom(biom_json_ids([1, 2], [10, 20]))
    assert list(tdata.obs_names) == ["1", "2"] and list(tdata.var_names) == ["10", "20"]


def test_read_biom_duplicate_ids_are_rejected_by_biom(biom_json_ids):
    # biom-format validates ids itself (biom/err.py SAMPDUP, default state "raise").
    with pytest.raises(TableException, match="Duplicate sample IDs"):
        bt.io.read_biom(biom_json_ids(["S1", "S1"], ["a", "b"]))


def test_read_biom_attaches_tree(biom_hdf5, newick):
    tree = bt.io.read_biom(biom_hdf5, tree=newick).vart["phylo"]
    assert {n for n in tree.nodes if tree.out_degree(n) == 0} == {"OTU_1", "OTU_2", "OTU_3", "OTU_4"}


def test_read_biom_tree_mismatch_keeps_shared_features(biom_hdf5, tmp_path):
    path = tmp_path / "partial.nwk"
    path.write_text("((OTU_1:1,OTU_2:1):1,OTU_9:1);")
    with pytest.warns(UserWarning, match="2 feature"):
        tdata = bt.io.read_biom(biom_hdf5, tree=path)
    assert list(tdata.var_names) == ["OTU_1", "OTU_2"]


def test_read_biom_without_taxonomy_has_no_rank_columns(tmp_path, biom_table):
    path = tmp_path / "bare.biom"
    bare = biom.Table(biom_table.matrix_data, biom_table.ids("observation"), biom_table.ids())
    with biom_open(str(path), "w") as handle:
        bare.to_hdf5(handle, "biotapy tests")
    assert not set(bt.io.read_biom(path).var.columns) & {"kingdom", "phylum"}


def test_read_biom_output_saves_to_h5td(biom_hdf5, tmp_path):
    tdata = bt.io.read_biom(biom_hdf5)
    assert tdata.var["species"].isna().all()
    expected = tdata.var.copy()  # the writer turns str columns into categoricals in place
    tdata.write_h5td(tmp_path / "x.h5td")
    back = treedata.read_h5td(tmp_path / "x.h5td")
    pd.testing.assert_frame_equal(back.var.astype("str"), expected)


def test_read_biom_records_counts_and_provenance(biom_hdf5):
    meta = bt.io.read_biom(biom_hdf5).uns["biotapy"]
    assert meta["x_kind"] == "counts" and '"io.read_biom"' in meta["provenance"][-1]


def test_read_biom_infers_relative_abundance(tmp_path):
    path = tmp_path / "relative.biom"
    table = biom.Table(np.array([[0.25, 0.5], [0.75, 0.5]]), ["OTU_1", "OTU_2"], ["S1", "S2"])
    with biom_open(str(path), "w") as handle:
        table.to_hdf5(handle, "biotapy tests")
    assert bt.io.read_biom(path).uns["biotapy"]["x_kind"] == "relative"


@pytest.mark.parametrize("fmt", ["hdf5", "json"])
def test_write_biom_round_trips_toy(tmp_path, fmt):
    toy = bt.datasets.toy()
    path = tmp_path / "toy.biom"
    bt.io.write_biom(toy, path, fmt=fmt)
    back = bt.io.read_biom(path)
    assert list(back.obs_names) == list(toy.obs_names) and list(back.var_names) == list(toy.var_names)
    np.testing.assert_array_equal(back.X.toarray(), toy.X.toarray())
    pd.testing.assert_frame_equal(back.var[RANK_COLUMNS].fillna("-"), toy.var[RANK_COLUMNS].fillna("-"))
    assert back.obs["group"].tolist() == toy.obs["group"].astype(str).tolist()


def test_write_biom_stores_features_by_samples(tmp_path):
    path = tmp_path / "toy.biom"
    bt.io.write_biom(bt.datasets.toy(), path)
    table = biom.load_table(str(path))
    assert table.shape == (8, 6) and list(table.ids("observation"))[:2] == ["f1", "f2"]


def test_write_biom_prefixes_ranks_in_canonical_order(tmp_path):
    path = tmp_path / "toy.biom"
    bt.io.write_biom(bt.datasets.toy(), path)
    table = biom.load_table(str(path))
    assert list(table.metadata("f1", "observation")["taxonomy"]) == [
        "k__Bacteria",
        "p__Firmicutes",
        "c__Clostridia",
        "o__Lachnospirales",
        "f__Lachnospiraceae",
        "g__Blautia",
    ]
    assert list(table.metadata("f8", "observation")["taxonomy"])[-1] == "g__"


def test_write_biom_without_taxonomy_writes_no_metadata(tmp_path):
    adata = ad.AnnData(
        X=sp.csr_matrix(np.eye(2)), obs=pd.DataFrame(index=["s1", "s2"]), var=pd.DataFrame(index=["a", "b"])
    )
    path = tmp_path / "bare.biom"
    bt.io.write_biom(adata, path)
    assert biom.load_table(str(path)).metadata(axis="observation") is None


def test_write_biom_leaves_input_alone(tmp_path, assert_unchanged):
    toy = bt.datasets.toy()
    before = toy.copy()
    bt.io.write_biom(toy, tmp_path / "toy.biom")
    assert_unchanged(before, toy)


def test_write_biom_does_not_write_the_tree(tmp_path):
    path = tmp_path / "toy.biom"
    bt.io.write_biom(bt.datasets.toy(), path)
    assert "phylo" not in bt.io.read_biom(path).vart


@pytest.mark.parametrize("fmt", ["hdf5", "json"])
def test_write_biom_round_trips_missing_metadata_as_nan(tmp_path, fmt):
    toy = bt.datasets.toy()
    toy.obs["site"] = ["gut", None, "gut", "skin", None, "skin"]
    path = tmp_path / "toy.biom"
    bt.io.write_biom(toy, path, fmt=fmt)
    back = bt.io.read_biom(path).obs["site"]
    assert back.isna().tolist() == [False, True, False, False, True, False]
    assert back.dropna().tolist() == ["gut", "gut", "skin", "skin"]


def _bare(counts, samples, features):
    return ad.AnnData(X=sp.csr_matrix(counts), obs=pd.DataFrame(index=samples), var=pd.DataFrame(index=features))


@pytest.mark.parametrize("fmt", ["hdf5", "json"])
def test_write_biom_round_trips_all_zero_sample_and_feature(tmp_path, fmt):
    counts = np.array([[3, 0, 1], [0, 0, 0], [2, 0, 5]])
    path = tmp_path / "zeros.biom"
    bt.io.write_biom(_bare(counts, ["s1", "s2", "s3"], ["a", "b", "c"]), path, fmt=fmt)
    back = bt.io.read_biom(path)
    assert list(back.obs_names) == ["s1", "s2", "s3"] and list(back.var_names) == ["a", "b", "c"]
    np.testing.assert_array_equal(back.X.toarray(), counts)


@pytest.mark.parametrize("fmt", ["hdf5", "json"])
def test_write_biom_round_trips_single_sample(tmp_path, fmt):
    path = tmp_path / "one.biom"
    bt.io.write_biom(_bare(np.array([[4, 0, 7]]), ["s1"], ["a", "b", "c"]), path, fmt=fmt)
    back = bt.io.read_biom(path)
    assert back.shape == (1, 3) and back.X.toarray().tolist() == [[4, 0, 7]]


_IDS = st.text(alphabet="abXY019_.", min_size=1, max_size=4)


@st.composite
def _tables(draw):
    n_samples, n_features = draw(st.integers(1, 4)), draw(st.integers(1, 6))
    counts = draw(
        st.lists(
            st.lists(st.integers(0, 20), min_size=n_features, max_size=n_features),
            min_size=n_samples,
            max_size=n_samples,
        )
    )
    depth = draw(st.integers(1, len(ALL_RANKS)))
    values = st.lists(st.sampled_from(["Alpha", "Beta_2", None]), min_size=n_features, max_size=n_features)
    features = draw(st.lists(_IDS, min_size=n_features, max_size=n_features, unique=True))
    var = pd.DataFrame({rank: draw(values) for rank in ALL_RANKS[:depth]}, index=features, dtype="str")
    samples = draw(st.lists(_IDS, min_size=n_samples, max_size=n_samples, unique=True))
    return ad.AnnData(X=sp.csr_matrix(np.array(counts)), obs=pd.DataFrame(index=samples), var=var)


@settings(max_examples=40, deadline=None)
@given(_tables())
def test_write_then_read_biom_preserves_counts_ids_and_ranks(adata):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "table.biom"
        bt.io.write_biom(adata, path)
        back = bt.io.read_biom(path)
    assert list(back.obs_names) == list(adata.obs_names) and list(back.var_names) == list(adata.var_names)
    np.testing.assert_array_equal(back.X.toarray(), adata.X.toarray())
    pd.testing.assert_frame_equal(back.var, adata.var)
