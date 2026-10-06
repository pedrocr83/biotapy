"""The rpy2 bridge for the da methods only R implements (decisions/r-bridge-before-ports); extra ``r``."""

import warnings
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


# Runs the function and returns its value, the messages of the R warnings it raised (muffled so that R does not print
# them again; message() output is not a warning condition and keeps rpy2's default handling) and the message of the
# error that stopped it. The caller's .Random.seed (or its absence) is put back on exit, so a method's set.seed
# leaves nothing in the user's R session (rules.md R3.4).
_CATCH_WARNINGS = """function(f) function(...) {
  seed <- globalenv()
  before <- get0(".Random.seed", seed, inherits = FALSE)
  on.exit(if (is.null(before)) suppressWarnings(rm(".Random.seed", envir = seed)) else assign(".Random.seed", before, seed))
  caught <- character()
  tryCatch({
    value <- withCallingHandlers(f(...), warning = function(cond) {
      caught <<- c(caught, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = caught, error = character())
  }, error = function(cond) list(value = NULL, warnings = caught, error = conditionMessage(cond)))
}"""


def r_function(code: str, *, package: str, func: str) -> Callable[..., Any]:
    """The R function ``code`` defines; ImportError naming the extra without rpy2, or the install line without ``package``.

    Calling it returns the R function's value, re-emits each R warning as a ``UserWarning`` and turns an R error into a
    ``RuntimeError``, both naming ``func``; the warnings raised before an error are listed in its message.
    It leaves R's random state as it found it.
    """
    robjects = import_optional("rpy2.robjects", extra="r")
    rpackages = import_optional("rpy2.robjects.packages", extra="r")
    try:
        # Loaded before the call: a package that `::` loads during it makes R print "stack imbalance" warnings.
        rpackages.importr(package)
    except rpackages.PackageNotInstalledError as err:
        msg = f'{func} needs the R package {package}. Install it in R with: BiocManager::install("{package}")'
        raise ImportError(msg) from err
    function = robjects.r(f"({_CATCH_WARNINGS})({code})")

    def run(*args: object) -> object:
        value, caught, error = function(*args).values()
        if len(error):
            before = f" (R warned first: {'; '.join(caught)})" if len(caught) else ""
            msg = f"{func}: R stopped: {error[0].strip()}{before}"
            raise RuntimeError(msg)
        for message in caught:
            warnings.warn(f"{func}: R warned: {message}", stacklevel=4)
        return value

    return run


def call_r(function: Callable[..., Any], /, *args: object) -> pd.DataFrame:
    """``function(*args)``: pandas arguments become R vectors and data frames, the R data frame it returns pandas."""
    robjects = import_optional("rpy2.robjects", extra="r")
    pandas2ri = import_optional("rpy2.robjects.pandas2ri", extra="r")
    # A local converter, never rpy2's global activation, so the user's own rpy2 session keeps its conversion rules.
    with (robjects.default_converter + pandas2ri.converter).context():
        # ModuleType attributes are Any, so the result needs a cast under mypy --strict (warn_return_any).
        return cast("pd.DataFrame", function(*args))
