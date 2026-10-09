---
type: Decision
title: Heavy dependencies are optional extras
description: torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily; `pip install biotapy` stays light.
tags: [packaging, dependencies]
status: stable
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T09:34:26Z }
commit: 0a34a3a
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Context
Target users range from notebook researchers to pharma pipelines. Pulling
torch or an R installation into every install is unacceptable.[^spec]

# Decision
- Core runtime deps, each added in the phase that first imports it:
  anndata, treedata (`>=0.3.1,<0.4`), networkx, numpy, scipy, pandas,
  scikit-bio (`>=0.7.4,<0.8`), pooch (Phase 0-1); scikit-learn
  (`>=1.8`, Phase 1, NMDS); matplotlib (`>=3.8`, Phase 1, `pl`); mudata
  (Phase 2). Phase 0 added numpy; the
  template adds `session-info2` (debug report referenced by the issue template).
  Phase 1 task 1.1 added scipy and pandas, plus the dev-only stubs
  `pandas-stubs` and `scipy-stubs` so `mypy --strict` can check them.
  Task 1.3 added treedata and networkx, plus the dev-only `types-networkx`.
  Task 1.7a added scikit-bio (`>=0.7.4,<0.8`), which brings biom-format,
  statsmodels and patsy. Task 1.7c declared biom-format (`>=2.1.16`) directly
  because `io` imports it. Task 1.10 added rdata (`>=1.1,<2`) and xarray, to
  read phyloseq objects saved from R natively
  ([phyloseq-import-route](phyloseq-import-route.md)). Task 1.11 added pooch,
  to download and cache `datasets.global_patterns`/`enterotype` from
  phyloseq's repository at run time instead of committing its (AGPL-3) data
  files (R6.6). Task 1.17 added scikit-learn (`>=1.8`, approved 2026-09-27):
  scikit-bio has no non-metric MDS, and 1.8 renamed `dissimilarity` to
  `metric`. It brings joblib, threadpoolctl and cloudpickle, and adds
  about 0.15 s to `import biotapy`. Task 1.18 added matplotlib (`>=3.8`,
  approved 2026-09-27) for `pl`. 3.8 is its first release with type
  information and CPython 3.12 wheels. `pl` imports it inside its functions,
  and pyplot only when it must make a figure, so `import biotapy` does not
  load it. It brings contourpy, cycler, fonttools, kiwisolver, pillow and
  pyparsing. Checkpoint D declared threadpoolctl (`>=3.5`, approved
  2026-10-03), already installed through scikit-learn: `tl.permanova` limits
  scikit-bio's OpenMP F-statistic to one thread with it. Phase 2 task 2.1b
  added mudata (`>=0.4`, approved 2026-10-03): a HUMAnN or PICRUSt2 table needs
  a community and a per-taxon modality over the same samples, which neither
  AnnData nor TreeData holds. It is pure Python (BSD-3); what it needs
  (`scverse-misc[settings]`, pydantic-settings, python-dotenv, pydantic) is
  already installed through anndata. Phase 3 task 3.6 added the extra `r`
  (`rpy2>=3.6.8`, approved 2026-10-05): ALDEx2 and MaAsLin 3 exist only in R.
  rpy2 (GPL-2.0-or-later) brings rpy2-rinterface, rpy2-robjects, cffi, jinja2
  and tzlocal. rpy2-rinterface has no Linux wheels on PyPI (`uv.lock`), so on
  Linux it builds from its sdist against the user's R and needs R's link
  headers (`libpcre2-dev`, `libdeflate-dev`, `libzstd-dev` among others on
  Ubuntu; without them rpy2 silently falls back to an ABI mode that cannot
  load R, so the `r-bridge` CI job sets `RPY2_CFFI_MODE=API` to fail the build
  instead, `.github/workflows/test.yaml`). biotapy imports it only through
  `import_optional` in `da/_r.py` and never bundles it. Phase 4 task 4.B0
  added the extra `torch` (`torch>=2.9`, approved 2026-10-09) for
  `ml.to_torch`: 2.9.0 is torch's first release with CPython 3.14 wheels.
  torch (BSD-3-Clause) brings filelock, fsspec, jinja2, networkx, setuptools,
  sympy (with mpmath) and typing-extensions; of these, fsspec, mpmath,
  setuptools and sympy are new to `uv.lock`, and no other version moves. uv
  installs it from PyTorch's CPU index (`[tool.uv]` in `pyproject.toml`:
  index `pytorch-cpu`, `https://download.pytorch.org/whl/cpu`, `explicit =
  true`, so no other package resolves there): `2.14.1+cpu` on Linux and
  Windows, `2.14.1` on macOS arm64 (no wheel exists for Intel macOS); the
  Linux wheel is 196 MB, against PyPI's 555 MB plus CUDA packages. The index
  is uv configuration only: the wheel's metadata says `torch>=2.9; extra ==
  'torch'`, and pip users get PyPI's torch.
- Extras (names fixed now so docs never change), each added in the phase that first uses it:

  | Extra | Pulls | First used |
  |---|---|---|
  | `numba` | numba (>=0.67, supports up to Python 3.14) | perf track, only on benchmark evidence; also unlocks scikit-bio's `engine="numba"` |
  | `r` | rpy2 (3.6.8) | Phase 3 `da` bridges |
  | `torch` | torch (>=2.9; under uv, CPU wheels from PyTorch's index) | Phase 4 `ml.to_torch` |
  | `plotnine` | plotnine | only if a `pl` function needs it |

- **`unifrac` (Striped UniFrac) is not an extra.** On 2026-09-26 PyPI had only
  an sdist that needs conda headers and cannot build on Windows; binaries are
  bioconda-only (linux-64, osx-64, osx-arm64). UniFrac therefore runs through
  scikit-bio. A conda-only `engine="striped"` may be added later, behind
  `import_optional`, if benchmarks justify it.

- Optional modules are imported inside the function through
  `biotapy._core.import_optional(name, extra)`, which raises `ImportError`
  naming the extra to install.
- A class that must inherit from an extra's base (rules.md R3.6), such as
  `ml.to_torch`'s torch `Dataset`, is defined inside the function after
  `import_optional` (`ml/_torch.py:to_torch`). mypy does not follow the
  extra (`follow_imports = "skip"`, `ignore_missing_imports` in
  `pyproject.toml`), so the type check gives the same answer whether or not
  the extra is installed.
- Adding any dependency (core or extra) requires the user's approval and a
  written reason in the PR.

# Rejected
- **Everything in core**: install size and conflict surface explode.
- **Separate plugin packages per extra**: premature; extras are enough until a
  plugin has its own release cadence.

# Consequences
- CI imports every module with no extras installed (`import-without-extras`),
  so a lazy import leaking to module level fails fast. A job with all extras
  comes with the first extra.
- `uv.lock` is not committed, so every CI job that runs `uv run` resolves
  afresh, and reads torch's versions from download.pytorch.org even when it
  installs no torch (measured: `uv lock` 0.10 s -> 0.37 s with a warm cache).

[^spec]: Python Microbiome Toolkit development report, section Module layout
