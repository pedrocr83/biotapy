import re
from collections.abc import Callable
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import as_csr

# Headless and identical to CI's MPLBACKEND=agg; tests/ is collected before the src/ doctests, so they get it too.
matplotlib.use("Agg")


def _assert_unchanged(before: AnnData, after: AnnData) -> None:
    assert (as_csr(before.X) != as_csr(after.X)).nnz == 0
    pd.testing.assert_frame_equal(before.obs, after.obs)
    pd.testing.assert_frame_equal(before.var, after.var)
    for slot in ("layers", "obsm", "obsp", "uns"):
        assert set(getattr(before, slot).keys()) == set(getattr(after, slot).keys()), slot
    for key in before.layers:
        assert (as_csr(before.layers[key]) != as_csr(after.layers[key])).nnz == 0, key
    before_meta = before.uns.get("biotapy", {})
    after_meta = after.uns.get("biotapy", {})
    # A stray write such as uns["biotapy"]["pcoa"] from a call with inplace=False.
    assert set(before_meta) == set(after_meta), "uns['biotapy']"
    assert before_meta.get("x_kind") == after_meta.get("x_kind")
    # provenance survives an h5 round-trip as an ndarray of str rather than a list.
    before_provenance = [str(entry) for entry in before_meta.get("provenance", [])]
    after_provenance = [str(entry) for entry in after_meta.get("provenance", [])]
    assert before_provenance == after_provenance


@pytest.fixture
def assert_unchanged() -> Callable[[AnnData, AnnData], None]:
    """Fail if a biotapy call mutated its input (rules.md R3.3)."""
    return _assert_unchanged


def _make_adata(dense: np.ndarray) -> AnnData:
    return AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )


# Session scope: Hypothesis rejects function-scoped fixtures in @given tests.
@pytest.fixture(scope="session")
def make_adata() -> Callable[[np.ndarray], AnnData]:
    """AnnData from a dense samples x features array, with samples ``s0..`` and features ``f0..``."""
    return _make_adata


# A MyST notebook's code cells and a page's fenced Python blocks, in page order.
_PAGE_CODE = re.compile(r"^```(?:\{code-cell\} ipython3|python)\n(.*?)^```$", re.MULTILINE | re.DOTALL)


def _run_page(path: Path) -> dict[str, object]:
    namespace: dict[str, object] = {}
    exec("\n".join(_PAGE_CODE.findall(path.read_text(encoding="utf-8"))), namespace)
    return namespace


@pytest.fixture(scope="session")
def run_page() -> Callable[[Path], dict[str, object]]:
    """Run a docs page's code, as a reader would, and return the names it defines (tests quoting a page's numbers)."""
    return _run_page
