"""Fixtures shared by `tests/` and the `src/biotapy` doctests (both testpaths in pyproject.toml).

Only conftest.py here is common to both trees; one under `src/biotapy/` would ship in the
wheel (rules.md R4), and this file is left out of both the wheel (hatchling packages only
`src/biotapy`) and the sdist (its explicit `build.targets.sdist.include` list, pyproject.toml).
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
