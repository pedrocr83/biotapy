import pytest

from biotapy._core import import_optional


def test_imports_an_installed_module():
    assert import_optional("json", extra="none").dumps([]) == "[]"


def test_missing_module_names_the_extra():
    with pytest.raises(ImportError, match=r"pip install 'biotapy\[torch\]'"):
        import_optional("biotapy_no_such_module", extra="torch")
