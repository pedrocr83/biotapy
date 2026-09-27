import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import treedata

import biotapy as bt

PHYLOSEQ = Path(__file__).parents[1] / "data" / "phyloseq"


def test_read_phyloseq_rds_matches_toy():
    tdata, toy = bt.io.read_phyloseq(PHYLOSEQ / "toy.rds"), bt.datasets.toy()
    assert list(tdata.obs_names) == list(toy.obs_names) and list(tdata.var_names) == list(toy.var_names)
    np.testing.assert_array_equal(tdata.X.toarray(), toy.X.toarray())
    assert tdata.obs["group"].astype(str).tolist() == toy.obs["group"].astype(str).tolist()
    assert tdata.var.loc["f1", "genus"] == "Blautia" and np.isnan(tdata.var.loc["f8", "genus"])
    tree = tdata.vart["phylo"]
    assert {n for n in tree.nodes if tree.out_degree(n) == 0} == {f"f{i}" for i in range(1, 9)}
    assert tdata.uns["biotapy"]["x_kind"] == "counts"


def test_read_phyloseq_single_object_rdata():
    assert bt.io.read_phyloseq(PHYLOSEQ / "toy.RData").shape == (6, 8)


def test_read_phyloseq_samples_as_rows_and_null_slots():
    tdata = bt.io.read_phyloseq(PHYLOSEQ / "samples_as_rows.rds")
    np.testing.assert_array_equal(tdata.X.toarray(), bt.datasets.toy().X.toarray())
    assert "phylo" not in tdata.vart and not set(tdata.var.columns) & {"kingdom", "genus"}
    assert "group" in tdata.obs.columns


def test_read_phyloseq_populated_refseq_raises_with_the_r_fix():
    with pytest.raises(ValueError, match=r"path=.*writeXStringSet"):
        bt.io.read_phyloseq(PHYLOSEQ / "with_refseq.rds")


def test_read_phyloseq_unparseable_file_names_the_cause_not_refseq(tmp_path):
    path = tmp_path / "not_r.RData"
    path.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError, match=r"path=.*cannot parse") as excinfo:
        bt.io.read_phyloseq(path)
    assert "refseq slot holds sequences" not in str(excinfo.value)


def test_read_phyloseq_two_objects_need_a_name():
    with pytest.raises(ValueError, match=r"name=.*toy.*toy_b"):
        bt.io.read_phyloseq(PHYLOSEQ / "two_objects.RData")
    assert bt.io.read_phyloseq(PHYLOSEQ / "two_objects.RData", name="toy_b").shape == (2, 8)


def test_read_phyloseq_unknown_name_is_named():
    with pytest.raises(KeyError, match="nope"):
        bt.io.read_phyloseq(PHYLOSEQ / "two_objects.RData", name="nope")


def test_read_phyloseq_rejects_a_non_phyloseq_file_and_names_what_it_holds():
    with pytest.raises(ValueError, match="path=.*phyloseq") as excinfo:
        bt.io.read_phyloseq(Path(__file__).parents[1] / "data" / "dada2" / "seqtab.rds")
    assert "DataArray" in str(excinfo.value)


def test_read_phyloseq_rds_content_with_rdata_suffix_reads_with_no_warning(tmp_path):
    path = tmp_path / "toy.RData"
    path.write_bytes((PHYLOSEQ / "toy.rds").read_bytes())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tdata = bt.io.read_phyloseq(path)
    assert not [w for w in caught if issubclass(w.category, UserWarning)]
    assert tdata.shape == (6, 8)


def test_read_phyloseq_rdata_content_with_rds_suffix_selects_by_name(tmp_path):
    path = tmp_path / "two.rds"
    path.write_bytes((PHYLOSEQ / "two_objects.RData").read_bytes())
    assert bt.io.read_phyloseq(path, name="toy_b").shape == (2, 8)


def test_read_phyloseq_name_on_a_single_object_file_is_named():
    with pytest.raises(ValueError, match="name="):
        bt.io.read_phyloseq(PHYLOSEQ / "toy.rds", name="toy")


def test_read_phyloseq_output_saves_to_h5td(tmp_path):
    tdata = bt.io.read_phyloseq(PHYLOSEQ / "toy.rds")
    assert tdata.obs["group"].dtype == pd.StringDtype(na_value=np.nan)
    expected_var = tdata.var.copy()
    tdata.write_h5td(tmp_path / "x.h5td")
    back = treedata.read_h5td(tmp_path / "x.h5td")
    pd.testing.assert_frame_equal(back.var.astype("str"), expected_var.astype("str"))
    assert back.obs["group"].astype(str).tolist() == tdata.obs["group"].astype(str).tolist()


def test_read_phyloseq_single_sample_keeps_all_zero_features():
    tdata = bt.io.read_phyloseq(PHYLOSEQ / "single_sample.rds")
    assert tdata.shape == (1, 8)
    np.testing.assert_array_equal(tdata.X.toarray(), bt.datasets.toy().X[:1].toarray())
    assert tdata.X.toarray()[0, [4, 7]].tolist() == [0, 0]
    tree = tdata.vart["phylo"]
    assert list(tdata.var_names) == [f"f{i}" for i in range(1, 9)]
    assert {n for n in tree.nodes if tree.out_degree(n) == 0} == {f"f{i}" for i in range(1, 9)}


def test_read_phyloseq_all_zero_sample_without_sample_data():
    tdata = bt.io.read_phyloseq(PHYLOSEQ / "zero_sample.rds")
    assert tdata.shape == (6, 8)
    assert tdata.X[5].nnz == 0
    np.testing.assert_array_equal(tdata.X.toarray()[:5], bt.datasets.toy().X.toarray()[:5])
    assert tdata.obs.columns.empty and list(tdata.obs_names) == [f"s{i}" for i in range(1, 7)]
    assert "phylo" not in tdata.vart
    assert tdata.uns["biotapy"]["x_kind"] == "counts"


def test_read_phyloseq_reads_uppercase_rds_extension(tmp_path):
    path = tmp_path / "toy.RDS"
    path.write_bytes((PHYLOSEQ / "toy.rds").read_bytes())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tdata = bt.io.read_phyloseq(path)
    assert not caught
    assert tdata.shape == (6, 8)
