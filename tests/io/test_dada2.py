import numpy as np
import pytest

import biotapy as bt

SEQTAB = '"","ACGTACGT","TTGACCAA","GGGCCCAA","CCCCAAAA"\n"S1",10,0,5,0\n"S2",0,0,0,0\n"S3",3,7,1,0\n'
TAXA = '"","Kingdom","Phylum","Genus"\n"ACGTACGT","Bacteria","Firmicutes","Blautia"\n"TTGACCAA","Bacteria","Bacteroidota",NA\n'


@pytest.fixture
def seqtab(tmp_path):
    path = tmp_path / "seqtab.csv"
    path.write_text(SEQTAB)
    return path


@pytest.fixture
def taxa(tmp_path):
    path = tmp_path / "taxa.csv"
    path.write_text(TAXA)
    return path


def test_read_dada2_names_asvs_and_keeps_sequences(seqtab):
    tdata = bt.io.read_dada2(seqtab)
    assert tdata.shape == (3, 4) and list(tdata.var_names) == ["ASV1", "ASV2", "ASV3", "ASV4"]
    assert tdata.var["sequence"].tolist() == ["ACGTACGT", "TTGACCAA", "GGGCCCAA", "CCCCAAAA"]
    np.testing.assert_array_equal(tdata.X.toarray()[0], [10, 0, 5, 0])


def test_read_dada2_keeps_all_zero_sample_and_feature(seqtab):
    tdata = bt.io.read_dada2(seqtab)
    assert tdata.X[1].nnz == 0 and tdata.X[:, 3].nnz == 0


def test_read_dada2_normalizes_taxa(seqtab, taxa):
    var = bt.io.read_dada2(seqtab, taxa).var
    assert {"kingdom", "phylum", "genus", "sequence"} <= set(var.columns)
    assert var.loc["ASV1", "genus"] == "Blautia" and np.isnan(var.loc["ASV2", "genus"])
    assert var.loc["ASV3", ["kingdom", "phylum", "genus"]].isna().all()


def test_read_dada2_reads_tsv(tmp_path):
    path = tmp_path / "seqtab.tsv"
    path.write_text(SEQTAB.replace(",", "\t"))
    assert bt.io.read_dada2(path).shape == (3, 4)


def test_read_dada2_single_sample(tmp_path):
    path = tmp_path / "one.csv"
    path.write_text('"","ACGTACGT"\n"S1",4\n')
    assert bt.io.read_dada2(path).shape == (1, 1)


def test_read_dada2_rejects_a_transposed_table(tmp_path):
    path = tmp_path / "asvs_by_samples.csv"
    path.write_text('"","S1","S2"\n"ACGTACGT",1,2\n')
    with pytest.raises(ValueError, match="samples x sequences"):
        bt.io.read_dada2(path)


def test_read_dada2_tree_tips_named_by_sequence(seqtab, tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("(((ACGTACGT:1,TTGACCAA:1):1,GGGCCCAA:1):1,CCCCAAAA:1);")
    phylo = bt.io.read_dada2(seqtab, tree=path).vart["phylo"]
    assert {n for n in phylo.nodes if phylo.out_degree(n) == 0} == {"ASV1", "ASV2", "ASV3", "ASV4"}
