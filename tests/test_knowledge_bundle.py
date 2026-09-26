"""OKF v0.2 conformance for .knowledge/ (rules.md R12.2)."""

from pathlib import Path

import pytest
import yaml

BUNDLE = Path(__file__).resolve().parents[1] / ".knowledge"
RESERVED = {"index.md", "log.md"}
CONCEPTS = sorted(p for p in BUNDLE.rglob("*.md") if p.name not in RESERVED)


def _frontmatter(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return {}
    block, found, _ = text[4:].partition("\n---\n")
    data = yaml.safe_load(block) if found else None
    return data if isinstance(data, dict) else {}


def _rel(path: Path) -> str:
    return path.relative_to(BUNDLE).as_posix()


def test_bundle_has_concepts():
    assert CONCEPTS


@pytest.mark.parametrize("path", CONCEPTS, ids=_rel)
def test_concept_has_a_type(path: Path):
    assert str(_frontmatter(path).get("type") or "").strip(), f"{_rel(path)} needs frontmatter with a non-empty `type`"


@pytest.mark.parametrize("path", CONCEPTS, ids=_rel)
def test_concept_is_listed_in_its_index(path: Path):
    index = path.parent / "index.md"
    assert index.is_file(), f"{_rel(path.parent)}/ has no index.md"
    assert f"({path.name})" in index.read_text(encoding="utf-8"), f"{_rel(path)} missing from {_rel(index)}"


def test_only_the_root_index_has_frontmatter():
    nested = [p for p in BUNDLE.rglob("index.md") if p.parent != BUNDLE]
    assert not [_rel(p) for p in nested if p.read_text(encoding="utf-8").startswith("---")]
