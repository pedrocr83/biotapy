from collections.abc import Callable

import pandas as pd
import pytest
from anndata import AnnData

from biotapy._core import as_csr


def _assert_unchanged(before: AnnData, after: AnnData) -> None:
    assert (as_csr(before.X) != as_csr(after.X)).nnz == 0
    pd.testing.assert_frame_equal(before.obs, after.obs)
    pd.testing.assert_frame_equal(before.var, after.var)
    for slot in ("layers", "obsm", "obsp", "uns"):
        assert set(getattr(before, slot).keys()) == set(getattr(after, slot).keys()), slot


@pytest.fixture
def assert_unchanged() -> Callable[[AnnData, AnnData], None]:
    """Fail if a biotapy call mutated its input (rules.md R3.3)."""
    return _assert_unchanged
