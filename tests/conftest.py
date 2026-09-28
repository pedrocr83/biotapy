from collections.abc import Callable

import matplotlib
import pandas as pd
import pytest
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
