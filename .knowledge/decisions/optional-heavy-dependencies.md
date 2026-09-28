---
type: Decision
title: Heavy dependencies are optional extras
description: torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily; `pip install biotapy` stays light.
tags: [packaging, dependencies]
status: stable
generated: { by: claude-code/claude-sonnet-5, at: 2026-09-28T07:47:57Z }
commit: 806bede
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
  pyparsing.
- Extras (names fixed now so docs never change), each added in the phase that first uses it:

  | Extra | Pulls | First used |
  |---|---|---|
  | `numba` | numba (>=0.67, supports up to Python 3.14) | perf track, only on benchmark evidence; also unlocks scikit-bio's `engine="numba"` |
  | `r` | rpy2 | Phase 3 `da` bridges |
  | `torch` | torch | Phase 4 `ml` loaders |
  | `plotnine` | plotnine | only if a `pl` function needs it |

- **`unifrac` (Striped UniFrac) is not an extra.** On 2026-09-26 PyPI had only
  an sdist that needs conda headers and cannot build on Windows; binaries are
  bioconda-only (linux-64, osx-64, osx-arm64). UniFrac therefore runs through
  scikit-bio. A conda-only `engine="striped"` may be added later, behind
  `import_optional`, if benchmarks justify it.

- Optional modules are imported inside the function through
  `biotapy._core.import_optional(name, extra)`, which raises `ImportError`
  naming the extra to install.
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

[^spec]: Python Microbiome Toolkit development report, section Module layout
