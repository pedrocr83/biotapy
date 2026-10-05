"""The generated Coming-from-R table covers phyloseq's 31 core functions (Phase 1 exit gate)."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("coming_from_r", ROOT / "docs" / "extensions" / "coming_from_r.py")
coming_from_r = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coming_from_r)

# .knowledge/roadmap/phase-1-core.md, "The 31 phyloseq functions the table must cover".
PHYLOSEQ_31 = """otu_table sample_data tax_table phy_tree refseq nsamples ntaxa sample_names taxa_names
sample_sums taxa_sums rank_names sample_variables get_taxa_unique prune_taxa prune_samples subset_taxa
subset_samples filter_taxa transform_sample_counts rarefy_even_depth tax_glom estimate_richness distance
UniFrac ordinate plot_bar plot_richness plot_ordination plot_heatmap import_biom""".split()
NOT_IN_0_2 = ["tip_glom", "merge_samples", "psmelt", "plot_tree", "plot_net"]


def test_the_list_has_31_functions():
    assert len(set(PHYLOSEQ_31)) == 31


@pytest.mark.parametrize("name", PHYLOSEQ_31)
def test_table_maps_each_of_the_31(name):
    cells = coming_from_r.rows()[f"phyloseq::{name}"]
    assert cells and "not in 0.2" not in cells


@pytest.mark.parametrize("name", NOT_IN_0_2)
def test_uncovered_functions_are_marked(name):
    assert coming_from_r.rows()[f"phyloseq::{name}"] == ["not in 0.2"]


@pytest.mark.parametrize(
    ("r_name", "function"), [("importHUMAnN", "read_humann"), ("importMetaPhlAn", "read_metaphlan")]
)
def test_mia_importers_map_to_the_function_readers(r_name, function):
    assert coming_from_r.rows()[f"mia::{r_name}"] == [f"{{func}}`bt.io.{function} <biotapy.io.{function}>`"]


def test_plot_functions_link_to_pl():
    assert coming_from_r.rows()["phyloseq::plot_bar"] == ["{func}`bt.pl.bar <biotapy.pl.bar>`"]


def test_render_is_one_sorted_markdown_table():
    lines = coming_from_r.render().splitlines()
    assert lines[:2] == ["| R | biotapy |", "|---|---|"]
    names = [line.split("`")[1] for line in lines[2:]]
    assert names == sorted(names, key=str.lower) and len(names) == len(coming_from_r.rows())


def test_an_idiom_a_docstring_already_maps_raises(tmp_path, monkeypatch):
    idioms = tmp_path / "r_idioms.toml"
    idioms.write_text('[idioms]\n"phyloseq::tax_glom" = "`x`"\n', encoding="utf-8")
    monkeypatch.setattr(coming_from_r, "IDIOMS", idioms)
    with pytest.raises(ValueError, match=r"\['phyloseq::tax_glom'\] are mapped by a docstring"):
        coming_from_r.rows()


def test_write_table_writes_the_file_the_page_includes(tmp_path):
    coming_from_r.write_table(SimpleNamespace(srcdir=str(tmp_path)))
    assert (tmp_path / coming_from_r.TABLE).read_text(encoding="utf-8") == coming_from_r.render()
    page = (ROOT / "docs" / "coming_from_r.md").read_text(encoding="utf-8")
    assert f"```{{include}} {coming_from_r.TABLE.as_posix()}\n```" in page
