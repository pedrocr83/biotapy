"""Fixtures shared by `tests/` and the `src/biotapy` doctests (both testpaths in pyproject.toml).

Only conftest.py here is common to both trees; one under `src/biotapy/` would ship in the
wheel (rules.md R4), and this file is left out of the wheel (hatchling packages only
`src/biotapy`) and shipped in the sdist (`build.targets.sdist.include`, pyproject.toml), whose
tests need its marker hook for the `to_torch` doctest.
"""

from collections.abc import Iterator

import matplotlib
import pytest

# Headless, and needed here too: a run that collects only src/biotapy (skipping tests/,
# whose conftest.py also sets this) would otherwise pick a display backend. A matching
# backend makes matplotlib.use a no-op, so setting it twice is harmless.
matplotlib.use("Agg")


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    """Close every figure a test or doctest left open (the pl doctests draw with ax=None)."""
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Give the doctests of a module that needs the extra `torch` its marker, before `-m` deselects."""
    for item in items:
        if isinstance(item, pytest.DoctestItem) and item.name.startswith("biotapy.ml._torch."):
            item.add_marker(pytest.mark.torch)
