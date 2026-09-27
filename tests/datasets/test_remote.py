from pathlib import Path

import pytest

import biotapy as bt
from biotapy.datasets import _remote

TOY_RDS = Path(__file__).parents[1] / "data" / "phyloseq" / "toy.rds"


def test_loaders_read_the_fetched_file(monkeypatch):
    # Monkeypatching a private helper is normally R4.9-only for `_core`; this is the
    # one exception (R11.4): it is the only way to test the loaders offline without
    # committing phyloseq's (AGPL-3) data files.
    fetched = []
    monkeypatch.setattr(_remote, "_fetch", lambda name: fetched.append(name) or str(TOY_RDS))
    assert bt.datasets.global_patterns().shape == (6, 8)
    assert bt.datasets.enterotype().shape == (6, 8)
    assert fetched == ["GlobalPatterns.RData", "enterotype.RData"]


@pytest.mark.network
def test_global_patterns_downloads_and_loads():
    tdata = bt.datasets.global_patterns()
    tree = tdata.vart["phylo"]
    assert tdata.shape == (26, 19216) and sum(1 for n in tree.nodes if tree.out_degree(n) == 0) == 19216
    assert tdata.uns["biotapy"]["x_kind"] == "counts" and "phylum" in tdata.var.columns


@pytest.mark.network
def test_enterotype_downloads_as_relative_abundance():
    tdata = bt.datasets.enterotype()
    assert tdata.shape == (280, 553) and tdata.uns["biotapy"]["x_kind"] == "relative"
    assert "phylo" not in tdata.vart
