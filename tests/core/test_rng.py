import numpy as np
import pytest

from biotapy._core import as_generator


def test_int_seed_is_reproducible():
    assert as_generator(7).integers(10**9) == as_generator(7).integers(10**9)


def test_generator_is_passed_through():
    rng = np.random.default_rng(1)
    assert as_generator(rng) is rng


def test_none_gives_a_generator():
    assert isinstance(as_generator(None), np.random.Generator)


def test_legacy_randomstate_is_rejected():
    with pytest.raises(TypeError, match="seed"):
        as_generator(np.random.RandomState(0))  # type: ignore[arg-type]
