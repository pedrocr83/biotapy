"""The rpy2 bridge for the da methods only R implements (decisions/r-bridge-before-ports); extra ``r``."""

from collections.abc import Callable
from typing import Any, cast

import numpy as np
import pandas as pd

from biotapy._core import as_generator, import_optional

# R's set.seed takes a 32-bit signed integer; drawing below 2**31 - 1 keeps it non-negative and never R's NA.
_SEED_BOUND = 2**31 - 1


def r_seed(seed: int | np.random.Generator | None) -> int:
    """The integer for R's ``set.seed``, drawn once from ``as_generator(seed)`` (rules.md R3.4)."""
    return int(as_generator(seed).integers(_SEED_BOUND))


def r_function(code: str, *, package: str, func: str) -> Callable[..., Any]:
    """The R function ``code`` defines; ImportError naming the extra without rpy2, or the install line without ``package``."""
    robjects = import_optional("rpy2.robjects", extra="r")
    rpackages = import_optional("rpy2.robjects.packages", extra="r")
    try:
        # Loaded before the call: a package that `::` loads during it makes R print "stack imbalance" warnings.
        rpackages.importr(package)
    except rpackages.PackageNotInstalledError as err:
        msg = f'{func} needs the R package {package}. Install it in R with: BiocManager::install("{package}")'
        raise ImportError(msg) from err
    return cast("Callable[..., Any]", robjects.r(code))


def call_r(function: Callable[..., Any], /, *args: object) -> pd.DataFrame:
    """``function(*args)``: pandas arguments become R vectors and data frames, the R data frame it returns pandas."""
    robjects = import_optional("rpy2.robjects", extra="r")
    pandas2ri = import_optional("rpy2.robjects.pandas2ri", extra="r")
    # A local converter, never rpy2's global activation, so the user's own rpy2 session keeps its conversion rules.
    with (robjects.default_converter + pandas2ri.converter).context():
        return cast("pd.DataFrame", function(*args))
