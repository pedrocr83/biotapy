from pathlib import Path

import pytest

TESTS = Path(__file__).parent
DATA_FILES = sorted(p for d in ("data", "golden") for p in (TESTS / d).rglob("*") if p.is_file())


@pytest.mark.parametrize("path", DATA_FILES, ids=lambda p: str(p.relative_to(TESTS)))
def test_data_files_stay_under_one_megabyte(path):
    # rules.md R6.6
    assert path.stat().st_size < 1_000_000


def test_versions_file_names_r_and_bioconductor():
    text = (TESTS / "golden" / "VERSIONS.txt").read_text()
    assert "R version 4.5.3" in text and "Bioconductor 3.22" in text
