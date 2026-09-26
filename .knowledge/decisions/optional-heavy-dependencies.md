---
type: Decision
title: Heavy dependencies are optional extras
description: torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily; `pip install biotapy` stays light.
tags: [packaging, dependencies]
status: stable
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T10:16:21Z }
commit: 3b29ffe
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
  scikit-bio (`>=0.7.4,<0.8`), matplotlib, pooch (Phase 0-1); scikit-learn
  (Phase 1, NMDS - pending approval); mudata (Phase 2). Phase 0 added numpy; the
  template adds `session-info2` (debug report referenced by the issue template).
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
- CI runs one job with all extras and one job with none, so a lazy import
  leaking to module level fails fast.

[^spec]: Python Microbiome Toolkit development report, section Module layout
