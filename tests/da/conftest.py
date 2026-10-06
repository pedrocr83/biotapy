import contextlib
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import biotapy as bt


# Session scope: GlobalPatterns takes seconds to load and glom; tests only read it.
@pytest.fixture(scope="session")
def benchmark():
    """The exit-gate data the golden files use: GlobalPatterns genera in >= 20% of samples, ``host`` and ``log_depth``."""
    tdata = bt.pp.filter_features(bt.pp.tax_glom(bt.datasets.global_patterns(), "genus"), min_prevalence=0.2)
    human = tdata.obs["SampleType"].isin(["Feces", "Skin", "Tongue"])
    tdata.obs["host"] = pd.Categorical(np.where(human, "human", "other"))
    tdata.obs["log_depth"] = np.log(np.asarray(tdata.X.sum(axis=1)).ravel())
    return tdata


class _Converter:
    """rpy2's converter as the bridge uses it: ``a + b``, then ``.context()``."""

    def __add__(self, other):
        return self

    def context(self):
        return contextlib.nullcontext()


@pytest.fixture
def fake_rpy2(monkeypatch):
    """The rpy2 modules the R bridges import, without R: every R function returns ``fake.output`` and logs its call.

    The bridges' own code (seeds, labels, orientation, the schema mapping) runs unchanged; only rpy2 is replaced.
    """
    fake = SimpleNamespace(output=None, installed={"ALDEx2", "maaslin3"}, calls=[], stderr=[])

    def r(code):
        def function(*args):
            fake.calls.append((code, args))
            for text in fake.stderr:
                callbacks.consolewrite_warnerror(text)
            if isinstance(fake.output, Exception):
                raise fake.output
            return fake.output

        return function

    callbacks = ModuleType("rpy2.rinterface_lib.callbacks")
    callbacks.consolewrite_warnerror = None
    embedded = ModuleType("rpy2.rinterface_lib.embedded")
    # rpy2's class for an R error; the bridge turns it into a RuntimeError.
    embedded.RRuntimeError = type("RRuntimeError", (RuntimeError,), {})
    robjects = ModuleType("rpy2.robjects")
    robjects.r = r
    robjects.default_converter = _Converter()
    packages = ModuleType("rpy2.robjects.packages")
    # rpy2's own class is an ImportError subclass too (rpy2.robjects.packages.LibraryError).
    packages.PackageNotInstalledError = type("PackageNotInstalledError", (ImportError,), {})

    def importr(name):
        if name not in fake.installed:
            raise packages.PackageNotInstalledError(f'The R package "{name}" is not installed.')

    packages.importr = importr
    pandas2ri = ModuleType("rpy2.robjects.pandas2ri")
    pandas2ri.converter = _Converter()
    modules = {"rpy2": ModuleType("rpy2"), "rpy2.robjects": robjects}
    modules |= {"rpy2.robjects.packages": packages, "rpy2.robjects.pandas2ri": pandas2ri}
    modules |= {"rpy2.rinterface_lib.callbacks": callbacks, "rpy2.rinterface_lib.embedded": embedded}
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)
    return fake
