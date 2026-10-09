"""The pseudocount step shared by pp's log-ratio transforms and ml's CLR."""

import numpy as np
import numpy.typing as npt

from ._matrix import as_csr, finite_non_negative
from ._warnings import warn_user


def check_pseudocount(pseudocount: object) -> None:
    """Raise unless ``pseudocount`` is a finite real number >= 0; a bool is refused."""
    if isinstance(pseudocount, bool) or not isinstance(pseudocount, int | float | np.integer | np.floating):
        msg = f"pseudocount must be a real number, got {pseudocount!r}"
        raise TypeError(msg)
    if not np.isfinite(pseudocount) or pseudocount < 0:
        msg = f"pseudocount must be a finite number >= 0, got {pseudocount!r}"
        raise ValueError(msg)


def pseudocounted(
    X: object, pseudocount: float, *, func: str, columns: npt.NDArray[np.intp] | None = None
) -> npt.NDArray[np.float64]:
    """``X`` (its ``columns``, in that order, if given) as a dense float64 array plus ``pseudocount``, checked > 0.

    ``X`` is anything :func:`as_csr` takes. ``func`` names the caller in messages.
    """
    check_pseudocount(pseudocount)
    matrix = as_csr(X).astype(np.float64)
    if columns is not None:
        # Reordered while sparse, so the dense copy below is the only one.
        matrix = matrix[:, columns]
    if not finite_non_negative(matrix):
        msg = f"{func} needs finite, non-negative values in X"
        raise ValueError(msg)
    positive = matrix.data[matrix.data > 0]
    if positive.size and pseudocount > positive.min():
        warn_user(
            f"pseudocount={pseudocount} is larger than the smallest non-zero value in X ({positive.min():.3g}), "
            f"so it swamps the rarest features; for relative abundances pass a pseudocount on their scale ({func})"
        )
    # scikit-bio's log-ratio functions need dense input (rules.md R6.2): one dense copy of X.
    values = matrix.toarray()
    values += pseudocount
    if np.any(values <= 0):
        msg = f"X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to {func}"
        raise ValueError(msg)
    return values
