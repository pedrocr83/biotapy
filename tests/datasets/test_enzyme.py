from pathlib import Path

import pytest

import biotapy as bt
from biotapy.datasets import _enzyme

# Excerpts of ENZYME (CC BY 4.0); tests/data/enzyme/NOTICE.txt gives the source and changes.
DATA = Path(__file__).parents[1] / "data" / "enzyme"


@pytest.fixture
def offline(monkeypatch):
    # As in test_remote.py: the only offline route to the loader is its private _fetch (R11.4).
    monkeypatch.setattr(_enzyme, "_fetch", lambda name: str(DATA / name))


def test_lists_every_ancestor_of_an_ec_number(offline):
    edges = bt.datasets.enzyme()
    rows = edges[edges["child"] == "1.1.1.1"]
    assert rows[["parent", "level"]].values.tolist() == [
        ["1.-.-.-", "class"],
        ["1.1.-.-", "subclass"],
        ["1.1.1.-", "subsubclass"],
    ]


def test_internal_ids_reach_the_levels_above_them(offline):
    edges = bt.datasets.enzyme()
    assert edges[edges["child"] == "2.7.1.-"]["parent"].tolist() == ["2.-.-.-", "2.7.-.-"]
    assert edges[edges["child"] == "2.-.-.-"].empty


def test_parents_carry_enzclass_names(offline):
    edges = bt.datasets.enzyme()
    names = edges.drop_duplicates("parent").set_index("parent")["parent_name"]
    assert names["1.-.-.-"] == "Oxidoreductases"
    assert names["3.2.1.-"] == "Glycosidases, i.e. enzymes hydrolyzing O- and S-glycosyl compounds"


def test_deleted_entries_keep_their_place(offline):
    assert "1.1.1.74" in set(bt.datasets.enzyme()["child"])


@pytest.mark.parametrize("child", ["1.1.1.5", "1.1.1.n4"])
def test_transferred_and_preliminary_entries_sit_under_their_class(offline, child):
    edges = bt.datasets.enzyme()
    assert edges[edges["child"] == child]["parent"].tolist() == ["1.-.-.-", "1.1.-.-", "1.1.1.-"]


def test_attrs_name_the_release_and_licence(offline):
    attrs = bt.datasets.enzyme().attrs
    assert "02-Sep-2026" in attrs["source"] and attrs["license"] == "CC BY 4.0"


def test_a_file_without_a_release_line_raises(monkeypatch, tmp_path):
    (tmp_path / "enzyme.dat").write_text("ID   1.1.1.1\nDE   alcohol dehydrogenase.\n//\n", encoding="utf-8")
    (tmp_path / "enzclass.txt").write_text((DATA / "enzclass.txt").read_text(encoding="utf-8"), encoding="utf-8")
    monkeypatch.setattr(_enzyme, "_fetch", lambda name: str(tmp_path / name))
    with pytest.raises(ValueError, match="enzyme.dat"):
        bt.datasets.enzyme()


def _write_files(tmp_path, monkeypatch, *, dat, classes):
    (tmp_path / "enzyme.dat").write_text(dat, encoding="utf-8")
    (tmp_path / "enzclass.txt").write_text(classes, encoding="utf-8")
    monkeypatch.setattr(_enzyme, "_fetch", lambda name: str(tmp_path / name))


def test_mixed_releases_raise(monkeypatch, tmp_path):
    classes = (DATA / "enzclass.txt").read_text(encoding="utf-8").replace("02-Sep-2026", "01-Jan-2026")
    _write_files(tmp_path, monkeypatch, dat=(DATA / "enzyme.dat").read_text(encoding="utf-8"), classes=classes)
    with pytest.raises(ValueError, match=r"enzyme\.dat.*02-Sep-2026.*enzclass\.txt.*01-Jan-2026"):
        bt.datasets.enzyme()


def test_enzclass_without_class_lines_raises(monkeypatch, tmp_path):
    dat = (DATA / "enzyme.dat").read_text(encoding="utf-8")
    _write_files(tmp_path, monkeypatch, dat=dat, classes="Release:     02-Sep-2026\n")
    with pytest.raises(ValueError, match="enzclass.txt"):
        bt.datasets.enzyme()


@pytest.mark.network
def test_enzyme_downloads_and_parses():
    edges = bt.datasets.enzyme()
    assert edges["child"].nunique() > 8000 and set(edges["level"]) == {"class", "subclass", "subsubclass"}
    assert edges.attrs["source"].startswith("ENZYME release ")
