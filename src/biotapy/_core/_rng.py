"""Single entry point for randomness (rules.md R3.4)."""

import numpy as np


def as_generator(seed: int | np.random.Generator | None) -> np.random.Generator:
    """Return a NumPy Generator for ``seed`` without touching global state."""
    if isinstance(seed, np.random.Generator):
        return seed
    if seed is None or isinstance(seed, int | np.integer):
        return np.random.default_rng(seed)
    msg = f"seed must be an int, a numpy Generator or None, got {type(seed).__name__}"
    raise TypeError(msg)
