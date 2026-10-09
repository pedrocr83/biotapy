from pathlib import Path

import pooch

from biotapy._core import make_pooch


def test_the_cache_is_biotapy_data_dir_when_it_is_set(monkeypatch, tmp_path):
    monkeypatch.setenv("BIOTAPY_DATA_DIR", str(tmp_path))
    cache = make_pooch("https://example.org/data/", {"a.txt": None}, urls={"a.txt": "https://example.org/a.txt"})
    assert Path(cache.abspath) == tmp_path
    assert cache.registry == {"a.txt": None} and cache.get_url("a.txt") == "https://example.org/a.txt"


def test_without_biotapy_data_dir_the_cache_is_pooch_s_per_user_directory(monkeypatch):
    monkeypatch.delenv("BIOTAPY_DATA_DIR", raising=False)
    cache = make_pooch("https://example.org/data/", {})
    assert Path(cache.abspath) == Path(pooch.os_cache("biotapy"))
