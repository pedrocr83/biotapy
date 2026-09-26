import biom
import numpy as np
import pandas as pd
import pytest
from biom.exception import TableException
from biom.util import biom_open

import biotapy as bt


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


def test_read_biom_records_counts_and_provenance(biom_hdf5):
    meta = bt.io.read_biom(biom_hdf5).uns["biotapy"]
    assert meta["x_kind"] == "counts" and '"io.read_biom"' in meta["provenance"][-1]
