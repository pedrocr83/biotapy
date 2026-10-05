"""Every public function's docstring keeps the contracts/function-shape skeleton (rules.md R8.2)."""

import importlib.util
import re
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "coming_from_r", Path(__file__).resolve().parents[1] / "docs" / "extensions" / "coming_from_r.py"
)
coming_from_r = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coming_from_r)
PUBLIC = coming_from_r.public_functions()


def test_every_public_subpackage_is_covered():
    assert {name.split(".")[1] for name, _ in PUBLIC} == {"da", "datasets", "fn", "io", "pl", "pp", "tl"}


@pytest.mark.parametrize(("name", "function"), PUBLIC, ids=[name for name, _ in PUBLIC])
def test_docstring_has_the_contract_sections(name, function):
    doc = function.__doc__ or ""
    coming_from_r.r_equivalents(doc)  # raises unless there is exactly one parseable line
    assert re.search(r"^\s*Notes\n\s*-----\n\s*R equivalent: .+\n\s*Guide: :doc:`/guide/\w+`$", doc, re.MULTILINE), name
    assert re.search(r"^\s*Examples\n\s*--------\n\s*>>> ", doc, re.MULTILINE), name


@pytest.mark.parametrize(
    "doc",
    [
        "Notes\n-----\nGuide: :doc:`/guide/x`\n",
        "R equivalent: ``phyloseq::a``\nR equivalent: ``phyloseq::b``\n",
        "R equivalent: phyloseq::a\n",
        "R equivalent: ``phyloseq::a``; ``vegan::b``\n",
    ],
    ids=["missing", "twice", "no backticks", "semicolon"],
)
def test_malformed_r_lines_raise(doc):
    with pytest.raises(ValueError, match="exactly one line"):
        coming_from_r.r_equivalents(doc)


def test_r_line_items_and_none():
    assert coming_from_r.r_equivalents("    R equivalent: ``phyloseq::ordinate``, ``ape::pcoa``\n") == [
        "phyloseq::ordinate",
        "ape::pcoa",
    ]
    assert coming_from_r.r_equivalents("R equivalent: none\n") == []
