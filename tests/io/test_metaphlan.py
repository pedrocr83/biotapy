import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import treedata as td
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt
from biotapy._core import RANKS

# tests/data/metaphlan/NOTICE.txt: a real MetaPhlAn 4.0.6 profile from HUMAnN's test data (MIT).
DEMO = Path(__file__).parents[1] / "data" / "metaphlan" / "demo_metaphlan_bugs_list.tsv"
LINEAGE = "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae"
# A default MetaPhlAn 4.2 profile (synthetic): four "#" lines, the header, UNCLASSIFIED first.
PROFILE_4_2 = (
    "#mpa_vJan25_CHOCOPhlAnSGB_202503\n#metaphlan s1.fastq -o s1_profile.tsv\n#1000 reads processed\n"
    "#SampleID\tMetaphlan_Analysis\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
    "UNCLASSIFIED\t-1\t20.0\t\n"
    "k__Bacteria\t2\t80.0\t\n"
    f"{LINEAGE}\t2|976|200643|171549|815\t80.0\t\n"
    f"{LINEAGE}|g__Bacteroides\t2|976|200643|171549|815|816\t80.0\t\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus\t2|976|200643|171549|815|816|28116\t50.0\t\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871\t2|976|200643|171549|815|816|28116|\t50.0\t\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_SGB1\t2|976|200643|171549|815|816|\t30.0\ts__Bacteroides_x\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_SGB1|t__SGB1\t2|976|200643|171549|815|816||\t30.0\t\n"
)


def write(tmp_path, text, name="table.tsv"):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_reads_the_leaf_clades_of_a_metaphlan_4_profile():
    # The file skips o__Corynebacteriales; its parent c__Actinomycetia must still not count as a leaf.
    tdata = bt.io.read_metaphlan(DEMO)
    assert tdata.obs_names.tolist() == ["demo_metaphlan_bugs_list"]
    assert tdata.var_names.tolist() == ["SGB2091", "SGB1871", "SGB2301", "SGB1814"]
    np.testing.assert_allclose(tdata.X.toarray(), [[0.5285687, 0.346816, 0.1214543, 0.0031609]], rtol=1e-12)
    assert tdata.uns["biotapy"]["x_kind"] == "relative"


def test_var_holds_the_ranks_of_each_leaf():
    var = bt.io.read_metaphlan(DEMO).var
    assert var.columns.tolist() == list(RANKS)
    assert var.loc["SGB1871"].tolist() == [
        "Bacteria",
        "Bacteroidetes",
        "Bacteroidia",
        "Bacteroidales",
        "Bacteroidaceae",
        "Bacteroides",
        "Bacteroides_ovatus",
    ]
    assert var.loc["SGB2091", "species"] == "Corynebacterium_SGB2091"


@pytest.mark.parametrize("rank", RANKS)
def test_glom_to_each_rank_matches_metaphlans_own_rows(rank):
    # MetaPhlAn's parity check (contracts/r-golden-parity, statement 8): a clade's row is the sum of its leaves.
    rows = pd.read_csv(DEMO, sep="\t", skiprows=5, index_col=0)["relative_abundance"]
    last = rows.index.str.split("|").str[-1]
    at_rank = last.str.startswith(f"{rank[0]}__")
    expected = (rows[at_rank] / 100).set_axis(last[at_rank].str[3:])
    out = bt.pp.tax_glom(bt.io.read_metaphlan(DEMO), rank)
    # MetaPhlAn 4.0.6 printed no o__Corynebacteriales row, so compare the clades it printed.
    got = pd.Series(out.X.toarray()[0], index=out.var[rank].to_numpy())[expected.index]
    # MetaPhlAn prints each percentage rounded to 5 decimals, so a parent and the sum of its leaves differ by up to
    # (number of leaves) x 5e-6 % (1.1e-5 % measured for Bacteroidaceae).
    pd.testing.assert_series_equal(got.sort_index(), expected.sort_index(), check_names=False, rtol=0, atol=1e-6)


def test_unclassified_is_kept_with_every_rank_missing(tmp_path):
    tdata = bt.io.read_metaphlan(write(tmp_path, PROFILE_4_2, "s1_profile.tsv"))
    assert tdata.obs_names.tolist() == ["s1"]
    assert tdata.var_names.tolist() == ["UNCLASSIFIED", "SGB1871", "SGB1"]
    assert tdata.var.loc["UNCLASSIFIED"].isna().all()
    np.testing.assert_allclose(tdata.X.toarray(), [[0.2, 0.5, 0.3]])


def test_reads_a_merged_metaphlan_4_table(tmp_path):
    text = (
        "#mpa_vJan25_CHOCOPhlAnSGB_202503\nclade_name\ts1\ts2_profile\n"
        "UNCLASSIFIED\t20.0\t0.0\nk__Bacteria\t80.0\t100.0\n"
        f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus\t80.0\t100.0\n"
        f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871\t80.0\t100.0\n"
    )
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.obs_names.tolist() == ["s1", "s2"]
    np.testing.assert_allclose(tdata.X.toarray(), [[0.2, 0.8], [0.0, 1.0]])


def test_reads_a_merged_metaphlan_3_table_and_drops_its_taxids(tmp_path):
    # MetaPhlAn 3 leaves are species; a leaf may stop above species, which is then NaN.
    text = (
        "#mpa_v30_CHOCOPhlAn_201901\nclade_name\tNCBI_tax_id\tA\tB\n"
        "k__Bacteria\t2\t100.0\t100.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae\t2|1239|186801|186802|186803\t100.0\t100.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia\t2|1239|186801|186802|186803|841\t40.0\t0.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia\t2|1239|186801|186802|186803|572511\t60.0\t100.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Blautia_obeum\t2|1239|186801|186802|186803|572511|40520\t60.0\t100.0\n"
    )
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.obs_names.tolist() == ["A", "B"]
    assert tdata.var_names.tolist() == ["Roseburia", "Blautia_obeum"]
    assert pd.isna(tdata.var.loc["Roseburia", "species"]) and tdata.var.loc["Roseburia", "genus"] == "Roseburia"


def test_reads_an_hmp2_style_table(tmp_path):
    # HMP2's taxonomic_profiles_3.tsv.gz: no "#" line, a "Feature\Sample" corner, UNKNOWN, CRLF line ends.
    text = (
        "Feature\\Sample\tCSM5FZ3N_P_profile\tCSM5FZ4M_profile\r\n"
        "UNKNOWN\t0\t100\r\nk__Bacteria\t100\t0\r\nk__Bacteria|p__Firmicutes\t100\t0\r\n"
    )
    tdata = bt.io.read_metaphlan(write(tmp_path, text, "taxonomic_profiles_3.tsv"))
    assert tdata.obs_names.tolist() == ["CSM5FZ3N_P", "CSM5FZ4M"]
    assert tdata.var_names.tolist() == ["UNKNOWN", "Firmicutes"]
    np.testing.assert_array_equal(tdata.X.toarray(), [[0.0, 1.0], [1.0, 0.0]])


def test_a_single_rank_table_reads_that_rank(tmp_path):
    # metaphlan --tax_lev s writes species rows only.
    text = "#mpa_v\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
    text += f"{LINEAGE}|g__B|s__B_a\t\t70.0\t\n{LINEAGE}|g__B|s__B_b\t\t30.0\t\n"
    assert bt.io.read_metaphlan(write(tmp_path, text)).var_names.tolist() == ["B_a", "B_b"]


def test_reads_read_stats_profiles_and_short_rows(tmp_path):
    # -t rel_ab_w_read_stats has a "-" coverage for UNCLASSIFIED; before 4.2.3 some rows lacked additional_species.
    text = (
        "#mpa_v\n#clade_name\tclade_taxid\trelative_abundance\tcoverage\testimated_number_of_reads_from_the_clade\n"
        "UNCLASSIFIED\t-1\t10.0\t-\t50\nk__Bacteria\t2\t90.0\n"
    )
    np.testing.assert_allclose(bt.io.read_metaphlan(write(tmp_path, text)).X.toarray(), [[0.1, 0.9]])


def test_reads_gzip(tmp_path):
    path = tmp_path / "demo_profile.tsv.gz"
    path.write_bytes(gzip.compress(DEMO.read_bytes()))
    tdata = bt.io.read_metaphlan(path)
    assert tdata.obs_names.tolist() == ["demo"] and tdata.n_vars == 4


def test_all_zero_sample_and_feature_are_kept(tmp_path):
    text = "clade_name\tA\tB\nk__Bacteria\t100.0\t0.0\nk__Bacteria|p__F\t100.0\t0.0\nk__Bacteria|p__G\t0.0\t0.0\n"
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.shape == (2, 2) and tdata.X[1].nnz == 0
    assert tdata.var_names.tolist() == ["F", "G"]


def test_an_empty_profile_has_no_features(tmp_path):
    text = "#mpa_v\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
    assert bt.io.read_metaphlan(write(tmp_path, text, "empty_profile.tsv")).shape == (1, 0)


@pytest.mark.parametrize(
    ("rows", "total"),
    [
        (f"k__Bacteria\t100.0\n{LINEAGE}|g__B\t60.0\n", 60.0),  # a leaf row removed
        ("d__Bacteria\t100.0\nd__Bacteria;p__Firmicutes\t100.0\n", 200.0),  # GTDB-style ";" lineages
    ],
    ids=["rows-removed", "semicolon-lineages"],
)
def test_leaves_that_do_not_sum_to_100_raise_naming_the_path(tmp_path, rows, total):
    path = write(tmp_path, "clade_name\tS1\n" + rows, "partial.tsv")
    with pytest.raises(ValueError, match=rf"partial\.tsv.*1 sample\(s\) do not sum to 100%: \{{'S1': {total}\}}"):
        bt.io.read_metaphlan(path)


@pytest.mark.parametrize(
    "text",
    [
        "clade_name\tS1\nk__Bacteria\tabc\n",
        "#mpa\n#clade_name\tNCBI_tax_id\trelative_abundance\nk__Bacteria\t2\t\n",
        "clade_name\tS1\nk__Bacteria\t50.0\t50.0\n",
        "clade_name\tS1\tS1\nk__Bacteria\t100.0\t100.0\n",
        "",
    ],
    ids=["not-a-number", "missing-value", "longer-row", "repeated-sample", "empty-file"],
)
def test_malformed_tables_raise_naming_the_path(tmp_path, text):
    with pytest.raises(ValueError, match=r"bad\.tsv"):
        bt.io.read_metaphlan(write(tmp_path, text, "bad.tsv"))


def test_repeated_leaf_names_raise_naming_the_path(tmp_path):
    text = "clade_name\tS1\nk__A|g__X\t50.0\nk__B|g__X\t50.0\n"
    with pytest.raises(ValueError, match=r"twice\.tsv.*duplicate var ids: \['X'\]"):
        bt.io.read_metaphlan(write(tmp_path, text, "twice.tsv"))


def test_round_trips_through_h5td(tmp_path):
    tdata = bt.io.read_metaphlan(DEMO)
    tdata.write_h5td(tmp_path / "m.h5td")
    back = td.read_h5td(tmp_path / "m.h5td")
    pd.testing.assert_frame_equal(back.var, tdata.var)


@given(st.lists(st.floats(0.01, 100, allow_nan=False), min_size=1, max_size=6))
def test_leaves_keep_every_percentage(weights):
    leaves = [round(100 * weight / sum(weights), 5) for weight in weights]
    clades = {f"{LINEAGE}|g__G{i % 2}|s__G{i % 2}_s{i}|t__SGB{i}": value for i, value in enumerate(leaves)}
    rows = dict(clades)
    for clade, value in clades.items():
        parts = clade.split("|")
        for depth in range(1, len(parts)):
            rows["|".join(parts[:depth])] = rows.get("|".join(parts[:depth]), 0) + value
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "t.tsv"
        path.write_text("clade_name\tS1\n" + "".join(f"{clade}\t{value!r}\n" for clade, value in rows.items()))
        tdata = bt.io.read_metaphlan(path)
    assert tdata.var_names.tolist() == [f"SGB{i}" for i in range(len(leaves))]
    np.testing.assert_allclose(tdata.X.toarray().ravel(), np.array(leaves) / 100, rtol=1e-12)


def test_a_sample_with_only_internal_rows_raises(tmp_path):
    # Sample B was profiled with --tax_lev g: its abundance sits on a row that is an ancestor in sample A.
    text = (
        "#mpa_vJan25\nclade_name\tA\tB\n"
        f"{LINEAGE}|g__Bacteroides\t100.0\t100.0\n"
        f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus\t100.0\t0.0\n"
        f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871\t100.0\t0.0\n"
    )
    with pytest.raises(ValueError, match=r"merged\.tsv.*1 sample\(s\) do not sum to 100%: \{'B': 0\.0\}"):
        bt.io.read_metaphlan(write(tmp_path, text, "merged.tsv"))


def test_var_names_carry_no_header_label(tmp_path):
    assert bt.io.read_metaphlan(DEMO).var.index.name is None
    text = "clade_name\tS1\nk__Bacteria\t100.0\n"
    assert bt.io.read_metaphlan(write(tmp_path, text)).var.index.name is None


def test_a_metaphlan_2_profile_is_named_after_its_file(tmp_path):
    text = "#SampleID\tMetaphlan2_Analysis\nk__Bacteria\t60.0\nk__Bacteria|p__Firmicutes\t60.0\nUNCLASSIFIED\t40.0\n"
    tdata = bt.io.read_metaphlan(write(tmp_path, text, "S9_profile.txt"))
    assert tdata.obs_names.tolist() == ["S9"]
    np.testing.assert_allclose(tdata.X.toarray(), [[0.6, 0.4]])


def test_a_merged_table_may_name_a_sample_relative_abundance(tmp_path):
    text = "#mpa_v\nclade_name\trelative_abundance\tS2\nk__Bacteria\t100.0\t100.0\n"
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.obs_names.tolist() == ["relative_abundance", "S2"]


def test_negative_abundances_raise_naming_the_path(tmp_path):
    text = "clade_name\tS1\nk__A|g__X\t150.0\nk__A|g__Y\t-50.0\n"
    with pytest.raises(ValueError, match=r"neg\.tsv.*negative"):
        bt.io.read_metaphlan(write(tmp_path, text, "neg.tsv"))


@pytest.mark.parametrize("kind", ["latin-1", "not-gzip"])
def test_unreadable_files_raise_naming_the_path(tmp_path, kind):
    path = tmp_path / ("bad.tsv.gz" if kind == "not-gzip" else "bad.tsv")
    path.write_bytes(b"clade_name\tS1\nk__Caf\xe9\t100.0\n")
    with pytest.raises(ValueError, match=r"bad\.tsv"):
        bt.io.read_metaphlan(path)


def test_read_humann_names_the_path_of_an_unreadable_file(tmp_path):
    path = tmp_path / "bad_genefamilies.tsv"
    path.write_bytes(b"# Gene Family\tS1_Abundance-RPKs\nK\xe91\t2.0\n")
    with pytest.raises(ValueError, match=r"bad_genefamilies\.tsv"):
        bt.io.read_humann(path)
