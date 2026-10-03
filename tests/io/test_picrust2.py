import gzip
import tempfile
from pathlib import Path

import mudata
import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt

# Synthetic PICRUSt2 2.6 outputs (PICRUSt2 is GPL-3: nothing is copied from it). Column names follow its
# documented headers. Abundances: S1 holds ASV1 (10) and 0042 (4), S2 holds 0042 (5.5) and RARE (2), S3 is empty.
# Copy numbers: ASV1 has EC:1.1.1.1 x1 and EC:2.7.1.1 x2; 0042 has EC:2.7.1.1 x1 and EC:3.2.1.1 x1; RARE has
# EC:1.1.1.1 x1. Each unstratified value is the sum of its contributions; EC:4.1.1.1 is predicted nowhere.
UNSTRAT = (
    "function\tS1\tS2\tS3\n"
    "EC:1.1.1.1\t10.0\t2.0\t0.0\n"
    "EC:2.7.1.1\t24.0\t5.5\t0.0\n"
    "EC:3.2.1.1\t4.0\t5.5\t0.0\n"
    "EC:4.1.1.1\t0.0\t0.0\t0.0\n"
)
CONTRIB_COLUMNS = (
    "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count"
    "\ttaxon_function_abun\ttaxon_rel_function_abun\tnorm_taxon_function_contrib\n"
)
CONTRIB = CONTRIB_COLUMNS + (
    "S1\tEC:1.1.1.1\tASV1\t10.0\t71.43\t1\t10.0\t71.43\t1.0\n"
    "S1\tEC:2.7.1.1\tASV1\t10.0\t71.43\t2\t20.0\t142.86\t0.833\n"
    "S1\tEC:2.7.1.1\t0042\t4.0\t28.57\t1\t4.0\t28.57\t0.167\n"
    "S1\tEC:3.2.1.1\t0042\t4.0\t28.57\t1\t4.0\t28.57\t1.0\n"
    "S2\tEC:1.1.1.1\tRARE\t2.0\t26.67\t1\t2.0\t26.67\t1.0\n"
    "S2\tEC:2.7.1.1\t0042\t5.5\t73.33\t1\t5.5\t73.33\t1.0\n"
    "S2\tEC:3.2.1.1\t0042\t5.5\t73.33\t1\t5.5\t73.33\t1.0\n"
)
# Pathways (8 columns): a pathway's contributions need not sum to its unstratified value.
PATH_UNSTRAT = "pathway\tS1\tS2\nPWY-1\t3.5\t1.0\nPWY-2\t0.5\t2.0\n"
PATH_CONTRIB = (
    "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count\ttaxon_function_abun\ttaxon_rel_function_abun\n"
    "S1\tPWY-1\tASV1\t10.0\t71.43\t0.5\t5.0\t35.71\n"
    "S1\tPWY-1\t0042\t4.0\t28.57\t0.25\t1.0\t7.14\n"
    "S2\tPWY-2\t0042\t5.5\t73.33\t0.5\t2.75\t36.67\n"
)


def write(tmp_path, text, name):
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.fixture
def mdata(tmp_path):
    unstrat = write(tmp_path, UNSTRAT, "pred_metagenome_unstrat.tsv")
    return bt.io.read_picrust2(unstrat, contrib=write(tmp_path, CONTRIB, "pred_metagenome_contrib.tsv"))


def test_reads_community_and_per_taxon_modalities(mdata):
    function, by_taxon = mdata["function"], mdata["function_by_taxon"]
    assert mdata.obs_names.tolist() == ["S1", "S2", "S3"]
    assert function.var_names.tolist() == ["1.1.1.1", "2.7.1.1", "3.2.1.1", "4.1.1.1"]
    assert by_taxon.var_names.tolist() == [
        "1.1.1.1|ASV1",
        "2.7.1.1|ASV1",
        "2.7.1.1|0042",
        "3.2.1.1|0042",
        "1.1.1.1|RARE",
    ]
    np.testing.assert_array_equal(by_taxon.X.toarray()[1], [0.0, 0.0, 5.5, 5.5, 2.0])


def test_var_matches_the_humann_layout(mdata):
    assert mdata["function"].var.columns.tolist() == ["name", "special"]
    var = mdata["function_by_taxon"].var
    assert var.columns.tolist() == ["function", "name", "taxon", "genus", "species", "special"]
    assert var["taxon"].tolist() == ["ASV1", "ASV1", "0042", "0042", "RARE"]
    assert var[["name", "genus", "species"]].isna().all().all() and not var["special"].any()


def test_x_kind_is_abundance_so_rarefy_refuses_it(mdata):
    assert mdata["function"].uns["biotapy"]["x_kind"] == "abundance"
    assert mdata["function_by_taxon"].uns["biotapy"]["x_kind"] == "abundance"
    with pytest.raises(ValueError, match="raw counts"):
        bt.pp.rarefy(mdata["function"], depth=5, seed=0)


def test_gene_family_contributions_sum_to_the_community_values(mdata):
    by_taxon = mdata["function_by_taxon"]
    summed = pd.DataFrame(by_taxon.X.toarray(), columns=by_taxon.var["function"]).T.groupby(level=0).sum().T
    community = pd.DataFrame(mdata["function"].X.toarray(), columns=mdata["function"].var_names)
    pd.testing.assert_frame_equal(summed, community[summed.columns], check_names=False)


def test_pathway_contributions_are_kept_as_written(tmp_path):
    unstrat = write(tmp_path, PATH_UNSTRAT, "path_abun_unstrat.tsv")
    mdata = bt.io.read_picrust2(unstrat, contrib=write(tmp_path, PATH_CONTRIB, "path_abun_contrib.tsv"))
    np.testing.assert_array_equal(mdata["function"].X.toarray(), [[3.5, 0.5], [1.0, 2.0]])
    np.testing.assert_array_equal(mdata["function_by_taxon"].X.toarray(), [[5.0, 1.0, 0.0], [0.0, 0.0, 2.75]])


def test_without_contrib_the_per_taxon_modality_is_empty(tmp_path):
    mdata = bt.io.read_picrust2(write(tmp_path, UNSTRAT, "unstrat.tsv"))
    assert mdata["function"].shape == (3, 4) and mdata["function_by_taxon"].shape == (3, 0)


def test_ec_ids_match_enzyme_style_hierarchies(mdata):
    edges = pd.DataFrame({"child": ["1.1.1.1", "2.7.1.1"], "parent": ["1.-.-.-", "2.-.-.-"], "level": "class"})
    out = bt.fn.func_glom(mdata["function"], "class", hierarchy=edges)
    assert out.var_names.tolist() == ["1.-.-.-", "2.-.-.-", "UNGROUPED"]


def test_ko_ids_are_kept(tmp_path):
    mdata = bt.io.read_picrust2(write(tmp_path, "function\tS1\nK00001\t1.5\n", "ko.tsv"))
    assert mdata["function"].var_names.tolist() == ["K00001"]


def test_all_zero_sample_and_feature_are_kept(mdata):
    assert mdata["function"].X[2].nnz == 0 and mdata["function_by_taxon"].X[2].nnz == 0
    assert mdata["function"].X[:, 3].nnz == 0


def test_single_sample(tmp_path):
    columns = CONTRIB_COLUMNS
    mdata = bt.io.read_picrust2(
        write(tmp_path, "function\tS1\nEC:1.1.1.1\t6.0\n", "u.tsv"),
        contrib=write(tmp_path, columns + "S1\tEC:1.1.1.1\tASV1\t3\t100\t2\t6.0\t200\t1\n", "c.tsv"),
    )
    assert mdata["function"].shape == (1, 1) and mdata["function_by_taxon"].shape == (1, 1)


def test_reads_gzip(tmp_path):
    unstrat, contrib = tmp_path / "u.tsv.gz", tmp_path / "c.tsv.gz"
    unstrat.write_bytes(gzip.compress(UNSTRAT.encode()))
    contrib.write_bytes(gzip.compress(CONTRIB.encode()))
    assert bt.io.read_picrust2(unstrat, contrib=contrib)["function_by_taxon"].n_vars == 5


def test_round_trips_through_h5mu(tmp_path, mdata):
    mdata.write_h5mu(tmp_path / "p.h5mu")
    back = mudata.read_h5mu(tmp_path / "p.h5mu")
    assert (back["function_by_taxon"].X != mdata["function_by_taxon"].X).nnz == 0
    assert back["function_by_taxon"].var["taxon"].tolist() == ["ASV1", "ASV1", "0042", "0042", "RARE"]


ROW = "\t1\t1\t1\t1.0\t1\t1\n"  # the six numeric contribution cells, taxon_function_abun = 1.0
# case -> (unstrat text, contrib text, the message expected)
MALFORMED = {
    "contrib-columns": (UNSTRAT, CONTRIB.replace("taxon_function_abun", "abun"), r"contrib=.*c\.tsv.*long-format"),
    "legacy-wide-contrib": (UNSTRAT, "function\tsequence\tS1\nEC:1.1.1.1\tASV1\t1.0\n", r"c\.tsv.*long-format"),
    "unknown-sample": (UNSTRAT, CONTRIB + "S9\tEC:1.1.1.1\tASV1" + ROW, r"c\.tsv.*samples that path lacks: \['S9'\]"),
    "unknown-function": (UNSTRAT, CONTRIB + "S1\tK00001\tASV1" + ROW, r"functions that path lacks: \['K00001'\]"),
    "repeated-row": (UNSTRAT, CONTRIB + "S1\tEC:1.1.1.1\tASV1" + ROW, r"c\.tsv.*repeats a sample, function and taxon"),
    "missing-taxon": (UNSTRAT, CONTRIB + "S1\tEC:1.1.1.1\t" + ROW, r"c\.tsv.*no function or taxon"),
    "non-number": (UNSTRAT, CONTRIB + "S1\tEC:4.1.1.1\tASV1\t1\t1\t1\tx\t1\t1\n", r"c\.tsv.*not a number"),
    "bar-in-taxon": (UNSTRAT, CONTRIB + "S1\tEC:4.1.1.1\tA|B" + ROW, r"path=.*u\.tsv.*contrib=.*c\.tsv.*one '\|'"),
    "repeated-sample": (UNSTRAT.replace("S3", "S2"), CONTRIB, r"path=.*u\.tsv.*repeats column names \['S2'\]"),
    "description": (
        "function\tdescription\tS1\nEC:1.1.1.1\tAlcohol dehydrogenase\t1.0\n",
        CONTRIB,
        r"u\.tsv.*'description'",
    ),
    "short-row": (UNSTRAT + "EC:1.1.1.1\t1\t1\n", CONTRIB, r"u\.tsv.*missing or NaN"),
}


@pytest.mark.parametrize("case", MALFORMED)
def test_malformed_inputs_raise_naming_the_file(tmp_path, case):
    unstrat, contrib, message = MALFORMED[case]
    with pytest.raises(ValueError, match=message):
        bt.io.read_picrust2(write(tmp_path, unstrat, "u.tsv"), contrib=write(tmp_path, contrib, "c.tsv"))


@given(
    st.lists(
        st.tuples(st.integers(0, 2), st.integers(0, 2), st.integers(0, 3), st.floats(0.01, 1e6, allow_nan=False)),
        min_size=1,
        max_size=12,
        unique_by=lambda row: row[:3],
    )
)
def test_each_contribution_lands_in_its_cell(rows):
    # (sample, function, taxon, value) rows in any order; each must appear once, at (sample, function|taxon).
    contrib = CONTRIB_COLUMNS + "".join(
        f"S{s}\tEC:{f}.1.1.1\tASV{t}\t1\t1\t1\t{value!r}\t1\t1\n" for s, f, t, value in rows
    )
    with tempfile.TemporaryDirectory() as tmp:
        unstrat = Path(tmp) / "u.tsv"
        unstrat.write_text("function\tS0\tS1\tS2\n" + "".join(f"EC:{f}.1.1.1\t0\t0\t0\n" for f in range(3)))
        (Path(tmp) / "c.tsv").write_text(contrib)
        by_taxon = bt.io.read_picrust2(unstrat, contrib=Path(tmp) / "c.tsv")["function_by_taxon"]
    X = pd.DataFrame(by_taxon.X.toarray(), index=by_taxon.obs_names, columns=by_taxon.var_names)
    got = [X.loc[f"S{s}", f"{f}.1.1.1|ASV{t}"] for s, f, t, _ in rows]
    # Looser than exact: pandas' default C float parser keeps about 15 significant digits (as in test_humann.py).
    np.testing.assert_allclose(got, [value for *_, value in rows], rtol=1e-12)
    assert by_taxon.X.nnz == len(rows)


@pytest.mark.parametrize("argument", ["path", "contrib"])
def test_a_file_that_is_not_utf8_text_raises_naming_its_argument(tmp_path, argument):
    bad = tmp_path / "bad.tsv"
    bad.write_bytes(b"function\tS1\n\xff\xfe\t1\n")
    unstrat, contrib = (bad, None) if argument == "path" else (write(tmp_path, UNSTRAT, "u.tsv"), bad)
    with pytest.raises(ValueError, match=rf"{argument}=.*bad\.tsv.*UTF-8"):
        bt.io.read_picrust2(unstrat, contrib=contrib)


@pytest.mark.parametrize("argument", ["path", "contrib"])
def test_an_empty_file_raises_naming_its_argument(tmp_path, argument):
    empty = write(tmp_path, "", "empty.tsv")
    unstrat, contrib = (empty, None) if argument == "path" else (write(tmp_path, UNSTRAT, "u.tsv"), empty)
    with pytest.raises(ValueError, match=rf"{argument}=.*empty\.tsv.*is empty"):
        bt.io.read_picrust2(unstrat, contrib=contrib)


def test_explicit_zero_contributions_are_not_stored(tmp_path):
    zero = CONTRIB + "S3\tEC:4.1.1.1\tASV1\t1\t1\t1\t0.0\t1\t1\n"
    by_taxon = bt.io.read_picrust2(write(tmp_path, UNSTRAT, "u.tsv"), contrib=write(tmp_path, zero, "c.tsv"))[
        "function_by_taxon"
    ]
    assert by_taxon.n_vars == 6 and by_taxon.X.nnz == 7


@pytest.mark.parametrize("argument", ["path", "contrib"])
def test_negative_values_raise_naming_their_argument(tmp_path, argument):
    unstrat = UNSTRAT.replace("10.0", "-10.0") if argument == "path" else UNSTRAT
    contrib = CONTRIB.replace("\t20.0\t142.86", "\t-20.0\t142.86") if argument == "contrib" else CONTRIB
    with pytest.raises(ValueError, match=rf"{argument}=.*negative"):
        bt.io.read_picrust2(write(tmp_path, unstrat, "u.tsv"), contrib=write(tmp_path, contrib, "c.tsv"))


def test_a_sample_with_a_nonzero_total_missing_from_contrib_raises(tmp_path):
    # PICRUSt2 drops zero rows only, so a sample with abundance always has contribution rows.
    only_s1 = "".join(line + "\n" for line in CONTRIB.splitlines()[:5])
    with pytest.raises(ValueError, match=r"contrib=.*c\.tsv.*no rows for samples.*\['S2'\]"):
        bt.io.read_picrust2(write(tmp_path, UNSTRAT, "u.tsv"), contrib=write(tmp_path, only_s1, "c.tsv"))


def test_an_all_zero_sample_may_be_missing_from_contrib(tmp_path):
    assert "S3" not in CONTRIB
    mdata = bt.io.read_picrust2(write(tmp_path, UNSTRAT, "u.tsv"), contrib=write(tmp_path, CONTRIB, "c.tsv"))
    assert mdata["function_by_taxon"].X[2].nnz == 0


def test_contrib_may_cover_a_subset_of_functions(tmp_path):
    # A pathway can be present at the community level with no contribution row for it.
    contrib = CONTRIB.replace("S1\tEC:3.2.1.1\t0042\t4.0\t28.57\t1\t4.0\t28.57\t1.0\n", "")
    contrib = contrib.replace("S2\tEC:3.2.1.1\t0042\t5.5\t73.33\t1\t5.5\t73.33\t1.0\n", "")
    mdata = bt.io.read_picrust2(write(tmp_path, UNSTRAT, "u.tsv"), contrib=write(tmp_path, contrib, "c.tsv"))
    assert "3.2.1.1" not in set(mdata["function_by_taxon"].var["function"])


def test_an_na_taxon_is_read_as_missing(tmp_path):
    contrib = CONTRIB + "S1\tEC:4.1.1.1\tNA\t1\t1\t1\t1.0\t1\t1\n"
    with pytest.raises(ValueError, match=r"c\.tsv.*no function or taxon"):
        bt.io.read_picrust2(write(tmp_path, UNSTRAT, "u.tsv"), contrib=write(tmp_path, contrib, "c.tsv"))


TRAITS = (
    "sequence\tEC:1.1.1.1\tEC:2.7.1.1\tEC:3.2.1.1\tmetadata_NSTI\n"
    "ASV1\t1\t2\t0\t0.03\n"
    "0042\t0\t1\t1\t0.12\n"
    "ASV9\t0\t0\t0\t1.5\n"
)


def test_traits_are_asvs_by_functions(tmp_path):
    traits = bt.io.read_picrust2_traits(write(tmp_path, TRAITS, "EC_predicted.tsv"))
    assert traits.index.tolist() == ["ASV1", "0042", "ASV9"]
    assert traits.columns.tolist() == ["1.1.1.1", "2.7.1.1", "3.2.1.1"]
    assert traits.dtypes.eq(np.float64).all() and traits.loc["ASV9"].eq(0).all()


def test_traits_read_gzip_and_keep_ko_ids(tmp_path):
    path = tmp_path / "KO_predicted.tsv.gz"
    path.write_bytes(gzip.compress(b"sequence\tK00001\tK00002\nASV1\t1\t0\n"))
    assert bt.io.read_picrust2_traits(path).columns.tolist() == ["K00001", "K00002"]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (TRAITS + "ASV1\t1\t1\t1\t0.1\n", r"repeats ASV ids: \['ASV1'\]"),
        (TRAITS + "ASV2\t1\tx\t1\t0.1\n", "not a number"),
        (TRAITS + "ASV2\t1\t1\n", "missing or NaN"),
        (TRAITS + "ASV2\t1\t-1\t1\t0.1\n", "negative"),
        ("", "is empty"),
    ],
    ids=["repeated-asv", "non-number", "short-row", "negative", "empty-file"],
)
def test_malformed_traits_raise_naming_the_path(tmp_path, text, message):
    with pytest.raises(ValueError, match=rf"traits\.tsv.*{message}"):
        bt.io.read_picrust2_traits(write(tmp_path, text, "traits.tsv"))


def test_an_extra_text_column_raises_naming_the_path(tmp_path):
    text = "sequence\tEC:1.1.1.1\tnote\nASV1\t1\thello\n"
    with pytest.raises(ValueError, match=r"traits\.tsv.*not a number"):
        bt.io.read_picrust2_traits(write(tmp_path, text, "traits.tsv"))


def test_functions_that_collide_after_the_ec_strip_raise_naming_the_path(tmp_path):
    text = "sequence\tEC:1.1.1.1\t1.1.1.1\nA\t1\t2\n"
    with pytest.raises(ValueError, match=r"traits\.tsv.*repeats function ids.*1\.1\.1\.1"):
        bt.io.read_picrust2_traits(write(tmp_path, text, "traits.tsv"))


@pytest.mark.parametrize("value", ["inf", "-inf"])
def test_non_finite_values_raise_naming_the_path(tmp_path, value):
    with pytest.raises(ValueError, match=r"traits\.tsv.*not finite"):
        bt.io.read_picrust2_traits(write(tmp_path, f"sequence\tEC:1.1.1.1\nA\t{value}\n", "traits.tsv"))
    with pytest.raises(ValueError, match=r"unstrat\.tsv.*not finite"):
        bt.io.read_picrust2(write(tmp_path, f"function\tS1\nEC:1.1.1.1\t{value}\n", "unstrat.tsv"))


def test_a_biom_converted_table_reads_its_hash_header(tmp_path):
    # biom convert --to-tsv (after qiime tools export) writes a "#" comment line, then a "#OTU ID" header.
    text = "# Constructed from biom file\n#OTU ID\tS1\tS2\nEC:1.1.1.1\t10.0\t2.0\nEC:2.7.1.1\t24.0\t5.5\nEC:3.2.1.1\t4.0\t5.5\n"
    mdata = bt.io.read_picrust2(write(tmp_path, text, "feature-table.tsv"))
    assert mdata.obs_names.tolist() == ["S1", "S2"]
    assert mdata["function"].var_names.tolist() == ["1.1.1.1", "2.7.1.1", "3.2.1.1"]
    np.testing.assert_array_equal(mdata["function"].X.toarray(), [[10.0, 24.0, 4.0], [2.0, 5.5, 5.5]])


def test_a_biom_converted_trait_table_raises_naming_its_first_cell(tmp_path):
    # q2-picrust2 exports only sample tables, so a trait table never has a "#OTU ID" header; before the shared
    # header rule this text was misread (functions '1' and '2'), with no error.
    text = "# Constructed from biom file\n#OTU ID\tEC:1.1.1.1\tEC:2.7.1.1\nASV1\t1\t2\nASV2\t0\t1\n"
    with pytest.raises(ValueError, match=r"traits\.tsv.*'#OTU ID'.*'sequence'.*read_picrust2\b"):
        bt.io.read_picrust2_traits(write(tmp_path, text, "traits.tsv"))


@pytest.mark.parametrize("first", ["function", "pathway", "#OTU ID", "OTU ID"])
def test_read_picrust2_accepts_picrust2_and_biom_headers(tmp_path, first):
    mdata = bt.io.read_picrust2(write(tmp_path, f"{first}\tS1\nEC:1.1.1.1\t1.0\n", "u.tsv"))
    assert mdata["function"].var_names.tolist() == ["1.1.1.1"]


@pytest.mark.parametrize(
    ("text", "cell"),
    [(TRAITS, "sequence"), ("\tS1\nEC:1.1.1.1\t1.0\n", ""), ("feature\tS1\nEC:1.1.1.1\t1.0\n", "feature")],
    ids=["trait-table", "empty-corner", "other"],
)
def test_read_picrust2_rejects_another_kind_of_table(tmp_path, text, cell):
    with pytest.raises(ValueError, match=rf"wrong\.tsv.*'{cell}'.*'function'.*read_picrust2_traits"):
        bt.io.read_picrust2(write(tmp_path, text, "wrong.tsv"))


@pytest.mark.parametrize(
    "text", [UNSTRAT, PATH_UNSTRAT, "\tEC:1.1.1.1\nASV1\t1\n"], ids=["ec", "pathway", "empty-corner"]
)
def test_read_picrust2_traits_rejects_another_kind_of_table(tmp_path, text):
    with pytest.raises(ValueError, match=r"wrong\.tsv.*'sequence'.*read_picrust2\b"):
        bt.io.read_picrust2_traits(write(tmp_path, text, "wrong.tsv"))
