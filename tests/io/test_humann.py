import gzip
import tempfile
from pathlib import Path

import mudata
import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt

# HUMAnN fixtures: tests/data/humann/NOTICE.txt says which are HUMAnN's (MIT) and which are synthetic.
DATA = Path(__file__).parents[1] / "data" / "humann"


def test_reads_community_and_stratified_rows():
    mdata = bt.io.read_humann(DATA / "genefamilies.tsv")
    function, by_taxon = mdata["function"], mdata["function_by_taxon"]
    assert function.var_names.tolist() == [
        "UNMAPPED",
        "UniRef90_A",
        "UniRef90_B",
        "UniRef90_C",
        "UniRef90_D",
        "UniRef90_unknown",
    ]
    assert by_taxon.n_vars == 6 and by_taxon.var["function"].tolist()[:2] == ["UniRef90_A", "UniRef90_A"]
    np.testing.assert_array_equal(function.X.toarray()[:, 1], [8.0, 2.0, 0.0])
    assert function.var.loc["UniRef90_A", "name"] == "alpha protein"


def test_samples_are_rows_with_the_unit_suffix_removed():
    mdata = bt.io.read_humann(DATA / "genefamilies.tsv")
    assert mdata["function"].obs_names.tolist() == ["S1", "S2", "S3"]
    assert mdata.obs_names.tolist() == ["S1", "S2", "S3"]


def test_specials_are_features_flagged_in_var():
    function = bt.io.read_humann(DATA / "pathabundance.tsv")["function"]
    assert function.var["special"].tolist() == [True, True, False, False]
    assert bt.io.read_humann(DATA / "pathabundance.tsv")["function_by_taxon"].var["special"].sum() == 2


@pytest.mark.parametrize(
    ("name", "kind"),
    [("genefamilies.tsv", "rpk"), ("gene_families.tsv", "cpm"), ("pathabundance.tsv", "abundance")],
)
def test_x_kind_comes_from_the_header(name, kind):
    mdata = bt.io.read_humann(DATA / name)
    assert mdata["function"].uns["biotapy"]["x_kind"] == kind
    assert mdata["function_by_taxon"].uns["biotapy"]["x_kind"] == kind


def test_renormalised_column_names_win_over_the_first_cell(tmp_path):
    path = tmp_path / "relab.tsv"
    path.write_text("# Gene Family HUMAnN v4.0.0.alpha.2 Adjusted CPMs\tS1-RELAB\nREADS_UNMAPPED\t0.25\nK1\t0.75\n")
    mdata = bt.io.read_humann(path)
    assert mdata["function"].uns["biotapy"]["x_kind"] == "relative"
    assert mdata.obs_names.tolist() == ["S1"]


def test_reads_humann_4_style_tables():
    # HUMAnN's own 4.x-style fixture: READS_UNMAPPED and bare strata such as "bug1".
    mdata = bt.io.read_humann(DATA / "gene_families.tsv")
    assert mdata["function"].var.loc["READS_UNMAPPED", "special"]
    assert mdata["function_by_taxon"].var["genus"].isna().all()
    assert mdata.obs_names.tolist() == ["HUMAnN_test"]


def test_reads_a_merged_humann_3_table():
    mdata = bt.io.read_humann(DATA / "multi_sample_genefamilies.tsv")
    assert mdata["function"].shape == (2, 25) and mdata["function_by_taxon"].shape == (2, 25)
    assert mdata["function_by_taxon"].var["species"].iloc[0] == "Dialister_invisus"


def test_community_rows_with_no_strata_give_an_empty_modality():
    mdata = bt.io.read_humann(DATA / "demo_pathabundance_with_names.tsv")
    assert mdata["function"].var_names.tolist() == ["UNMAPPED", "UNINTEGRATED"]
    assert mdata["function_by_taxon"].shape == (1, 0)


def test_reads_gzip(tmp_path):
    path = tmp_path / "genefamilies.tsv.gz"
    path.write_bytes(gzip.compress((DATA / "genefamilies.tsv").read_bytes()))
    assert bt.io.read_humann(path)["function"].n_vars == 6


def test_single_sample_table(tmp_path):
    path = tmp_path / "one.tsv"
    path.write_text("# Pathway\tS1_Abundance\nUNMAPPED\t1.5\nPWY-1: x\t2.5\nPWY-1: x|unclassified\t2.5\n")
    mdata = bt.io.read_humann(path)
    assert mdata["function"].shape == (1, 2) and mdata["function_by_taxon"].shape == (1, 1)


def test_all_zero_sample_and_feature_are_kept(tmp_path):
    path = tmp_path / "zeros.tsv"
    path.write_text("# Gene Family\tS1-RPKs\tS2-RPKs\nK1\t0.0\t3.0\nK2\t0.0\t0.0\nK2|unclassified\t0.0\t0.0\n")
    mdata = bt.io.read_humann(path)
    assert mdata["function"].shape == (2, 2) and mdata["function"].X[0].nnz == 0
    assert mdata["function_by_taxon"].var_names.tolist() == ["K2|unclassified"]


@given(st.lists(st.floats(0, 1e6, allow_nan=False, width=32), min_size=1, max_size=8))
def test_values_round_trip_through_a_written_table(values):
    rows = "".join(f"K{i}\t{value!r}\nK{i}|unclassified\t{value!r}\n" for i, value in enumerate(values))
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "t.tsv"
        path.write_text("# Gene Family\tS1_Abundance-RPKs\n" + rows)
        mdata = bt.io.read_humann(path)
    # Looser than exact: pandas' default C float parser keeps about 15 significant digits (5.1e-14
    # relative measured); exact parsing costs 2.7x on the 91 MB HMP2 table (2.1 s vs 5.8 s).
    np.testing.assert_allclose(mdata["function"].X.toarray().ravel(), values, rtol=1e-12)
    np.testing.assert_allclose(mdata["function_by_taxon"].X.toarray().ravel(), values, rtol=1e-12)


def test_round_trips_through_h5mu(tmp_path):
    mdata = bt.io.read_humann(DATA / "genefamilies.tsv")
    mdata.write_h5mu(tmp_path / "g.h5mu")
    back = mudata.read_h5mu(tmp_path / "g.h5mu")
    assert (back["function_by_taxon"].X != mdata["function_by_taxon"].X).nnz == 0


def test_without_a_hash_line_the_first_line_is_the_header(tmp_path):
    # The HMP2 merged tables: no "#", sample columns named after the joined files.
    path = tmp_path / "pathabundances_3.tsv"
    path.write_text("Feature\\Sample\tA_P_pathabundance_cpm\tB_P_pathabundance_cpm\nUNMAPPED\t334835\t314137\n")
    mdata = bt.io.read_humann(path)
    assert mdata.obs_names.tolist() == ["A_P", "B_P"]
    assert mdata["function"].uns["biotapy"]["x_kind"] == "cpm"


def test_pathway_coverage_raises(tmp_path):
    path = tmp_path / "cov.tsv"
    path.write_text("# Pathway\tS1_Coverage\nPWY-1\t0.5\n")
    with pytest.raises(ValueError, match="coverage"):
        bt.io.read_humann(path)


def test_repeated_samples_after_suffix_removal_raise(tmp_path):
    path = tmp_path / "dup.tsv"
    path.write_text("# Gene Family\tS1_Abundance-RPKs\tS1-RPKs\nK1\t1.0\t2.0\n")
    with pytest.raises(ValueError, match="duplicate obs ids"):
        bt.io.read_humann(path)


def test_header_only_table_has_no_features(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text("# Pathway\tS1_Abundance\tS2_Abundance\n")
    mdata = bt.io.read_humann(path)
    assert mdata["function"].shape == (2, 0) and mdata["function_by_taxon"].shape == (2, 0)


def test_empty_file_raises_naming_the_path(tmp_path):
    path = tmp_path / "empty_file.tsv"
    path.write_text("")
    with pytest.raises(ValueError, match=r"empty_file\.tsv.*is empty"):
        bt.io.read_humann(path)


def test_non_numeric_value_raises_naming_the_path(tmp_path):
    path = tmp_path / "text_value.tsv"
    path.write_text("# Gene Family\tS1-RPKs\nK1\tabc\nK2\t1.0\n")
    with pytest.raises(ValueError, match=r"text_value\.tsv.*not a number"):
        bt.io.read_humann(path)


@pytest.mark.parametrize(
    "rows",
    ["PWY\t3\t4\n", "K0\t1\nPWY\t3\t4\n", "K0\t1\nPWY\t\n", "K0\t1\nPWY\n", "K0\t1\t\n"],
    ids=["first-row-longer", "later-row-longer", "empty-cell", "row-shorter", "trailing-empty-cell"],
)
def test_ragged_or_missing_cells_raise_naming_the_path(tmp_path, rows):
    path = tmp_path / "ragged.tsv"
    path.write_text("# Pathway\tS1\n" + rows)
    with pytest.raises(ValueError, match=r"ragged\.tsv"):
        bt.io.read_humann(path)


@pytest.mark.parametrize(
    ("text", "error"),
    [
        ("# Gene Family\tS1_Abundance-RPKs\tS1-RPKs\nK1\t1.0\t2.0\n", "duplicate obs ids"),
        ("# Gene Family\tS1_Abundance-RPKs\nK1\t1.0\nK1\t2.0\n", "duplicate var ids"),
        ("# Gene Family\tS1_Abundance-RPKs\nK1|g__A.s__A_b|x\t1.0\n", r"one '\|'"),
    ],
    ids=["repeated samples", "repeated rows", "two bars"],
)
def test_id_errors_name_the_path(tmp_path, text, error):
    path = tmp_path / "ids.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match=rf"path='.*ids\.tsv': .*{error}") as info:
        bt.io.read_humann(path)
    assert isinstance(info.value.__cause__, ValueError)


def test_repeated_sample_columns_raise_naming_the_path(tmp_path):
    # pandas would rename the second "S1" to "S1.1" and read two samples.
    path = tmp_path / "twice.tsv"
    path.write_text("# Pathway\tS1\tS1\nPWY-1\t1.0\t2.0\n")
    with pytest.raises(ValueError, match=r"twice\.tsv.*repeats column names \['S1'\]"):
        bt.io.read_humann(path)


def test_row_without_an_id_raises_naming_the_path(tmp_path):
    path = tmp_path / "no_id.tsv"
    path.write_text("# Pathway\tS1\nPWY-1\t1.0\n\t2.0\n")
    with pytest.raises(ValueError, match=r"no_id\.tsv.*no id"):
        bt.io.read_humann(path)


@pytest.mark.parametrize("blank", [" ", "  "])
def test_row_with_a_blank_id_raises_naming_the_path(tmp_path, blank):
    path = tmp_path / "blank_id.tsv"
    path.write_text(f"# P\tS1\nA\t1\n{blank}\t2\n")
    with pytest.raises(ValueError, match=r"blank_id\.tsv.*no id"):
        bt.io.read_humann(path)


def test_header_with_an_empty_column_name_raises_naming_the_path(tmp_path):
    path = tmp_path / "trailing_tab.tsv"
    path.write_text("# P\tS1\t\nA\t1\t2\n")
    with pytest.raises(ValueError, match=r"trailing_tab\.tsv.*empty column name"):
        bt.io.read_humann(path)


def test_an_empty_corner_cell_is_read_as_an_unnamed_index(tmp_path):
    # pandas to_csv(sep="\t") and R write.table(col.names=NA) leave the header's first cell empty.
    path = tmp_path / "corner.tsv"
    path.write_text("\tS1_Abundance-RPKs\nUNMAPPED\t1.0\nK1|g__A.s__A_b\t2.0\n")
    mdata = bt.io.read_humann(path)
    assert mdata.obs_names.tolist() == ["S1"] and mdata["function"].uns["biotapy"]["x_kind"] == "rpk"
    assert mdata["function"].var_names.tolist() == ["UNMAPPED"]
    assert mdata["function_by_taxon"].var_names.tolist() == ["K1|g__A.s__A_b"]


def test_negative_abundance_raises_naming_the_path(tmp_path):
    path = tmp_path / "neg.tsv"
    path.write_text("# P\tS1\nA\t-1\n")
    with pytest.raises(ValueError, match=r"neg\.tsv.*negative"):
        bt.io.read_humann(path)
