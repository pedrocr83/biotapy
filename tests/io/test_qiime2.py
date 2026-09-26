import numpy as np
import pytest

import biotapy as bt

NEWICK = "((OTU_1:0.1,OTU_2:0.2):0.05,(OTU_3:0.3,OTU_4:0.4):0.1);"
TAXONOMY = (
    "Feature ID\tTaxon\tConfidence\n"
    "OTU_1\td__Bacteria; p__Firmicutes; c__Clostridia\t0.98\n"
    "OTU_2\td__Bacteria; p__Firmicutes\t0.7\n"
    "OTU_3\tUnassigned\t0.5\n"
)
METADATA = (
    "# written by hand\n"
    "sample-id\tdepth\tsite\n"
    "#q2:types\tnumeric\tcategorical\n"
    "S1\t1000\tgut\n"
    "S2\t\tskin\n"
    "\n"
    "S3\t500\tgut\n"
    "S9\t1\tgut\n"
)


@pytest.fixture
def table_qza(make_qza, biom_hdf5):
    return make_qza("table", "feature-table.biom", biom_hdf5.read_bytes())


def test_read_qiime2_table(table_qza):
    tdata = bt.io.read_qiime2(table_qza)
    assert tdata.shape == (3, 4) and list(tdata.obs_names) == ["S1", "S2", "S3"]
    assert tdata.uns["biotapy"]["x_kind"] == "counts"


def test_read_qiime2_taxonomy(table_qza, make_qza):
    taxonomy = make_qza("taxonomy", "taxonomy.tsv", TAXONOMY.encode())
    var = bt.io.read_qiime2(table_qza, taxonomy=taxonomy).var
    assert var.loc["OTU_1", ["kingdom", "phylum", "class"]].tolist() == ["Bacteria", "Firmicutes", "Clostridia"]
    assert var.loc["OTU_3", "kingdom"] == "Unassigned"
    assert var.loc["OTU_4"].isna().all()
    assert var["confidence"].tolist()[:3] == [0.98, 0.7, 0.5]


def test_read_qiime2_taxonomy_without_confidence(table_qza, make_qza):
    taxonomy = make_qza("taxonomy", "taxonomy.tsv", b"Feature ID\tTaxon\nOTU_1\tk__Bacteria\n")
    var = bt.io.read_qiime2(table_qza, taxonomy=taxonomy).var
    assert "confidence" not in var.columns and var.loc["OTU_1", "kingdom"] == "Bacteria"


def test_read_qiime2_tree(table_qza, make_qza):
    tree = make_qza("tree", "tree.nwk", NEWICK.encode())
    phylo = bt.io.read_qiime2(table_qza, tree=tree).vart["phylo"]
    assert {n for n in phylo.nodes if phylo.out_degree(n) == 0} == {"OTU_1", "OTU_2", "OTU_3", "OTU_4"}


def test_read_qiime2_wrong_artifact_names_the_argument(make_qza):
    taxonomy = make_qza("taxonomy", "taxonomy.tsv", TAXONOMY.encode())
    with pytest.raises(ValueError, match=r"table=.*feature-table\.biom"):
        bt.io.read_qiime2(taxonomy)


def test_read_qiime2_metadata(table_qza, tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text(METADATA)
    obs = bt.io.read_qiime2(table_qza, metadata=path).obs
    assert list(obs.index) == ["S1", "S2", "S3"]
    assert obs["site"].tolist() == ["gut", "skin", "gut"]
    assert obs["depth"].dtype.kind == "f" and np.isnan(obs.loc["S2", "depth"])


@pytest.mark.parametrize("header", ["id", "SampleID", "Sample-ID", "#SampleID", "sample_name"])
def test_read_qiime2_metadata_id_headers(table_qza, tmp_path, header):
    path = tmp_path / "metadata.tsv"
    path.write_text(f"{header}\tdepth\nS1\t10\nS2\t20\nS3\t30\n")
    assert bt.io.read_qiime2(table_qza, metadata=path).obs["depth"].tolist() == [10, 20, 30]


def test_read_qiime2_metadata_legacy_header_is_case_sensitive(table_qza, tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("#sampleid\tdepth\nS1\t10\n")
    with pytest.raises(ValueError, match="metadata="):
        bt.io.read_qiime2(table_qza, metadata=path)


def test_read_qiime2_metadata_row_longer_than_header_is_named(table_qza, tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("sample-id\tdepth\nS1\t10\textra\n")
    with pytest.raises(ValueError, match="metadata="):
        bt.io.read_qiime2(table_qza, metadata=path)


def test_read_qiime2_metadata_duplicate_ids_are_named(table_qza, tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("sample-id\tdepth\nS1\t10\nS1\t20\n")
    with pytest.raises(ValueError, match="S1"):
        bt.io.read_qiime2(table_qza, metadata=path)


def test_read_qiime2_metadata_ignores_empty_trailing_cells(table_qza, tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("sample-id\tdepth\nS1\t10\t\nS2\t20\t\nS3\t30\t\n")
    obs = bt.io.read_qiime2(table_qza, metadata=path).obs
    assert obs["depth"].tolist() == [10, 20, 30]


def test_read_qiime2_metadata_ignores_whitespace_trailing_cells(table_qza, tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("sample-id\tdepth\nS1\t10\t  \nS2\t20\t \t\nS3\t30\n")
    obs = bt.io.read_qiime2(table_qza, metadata=path).obs
    assert obs["depth"].tolist() == [10, 20, 30]
