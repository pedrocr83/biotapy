---
type: Contract
title: R golden parity
description: Every function with an R equivalent is tested against parquet golden files exported from pinned R; deterministic outputs match numerically, stochastic outputs match invariants.
tags: [testing, r, validation]
status: stable
paths: ["tests/r/**", "tests/golden/**", "tests/**/test_*.py"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: b77a226
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Statement
1. Golden files are produced by `tests/r/export_golden.R`, run only inside the
   pinned container `tests/r/Dockerfile`, never in normal CI.
2. Output: `tests/golden/<dataset>/<function>.parquet`, samples as rows
   (see [samples-as-rows](/decisions/samples-as-rows.md)), plus
   `tests/golden/VERSIONS.txt` listing R and package versions.
3. Regenerate only when the pinned R/Bioconductor version is bumped; the bump
   is its own commit.
4. Comparison rules by output kind:

   | Output | Compare | Default tolerance |
   |---|---|---|
   | Deterministic numeric (glom sums, relative, alpha, Bray-Curtis, UniFrac) | elementwise | `rtol=1e-7` |
   | PCoA coordinates | per axis, up to sign flip; eigenvalues elementwise | `rtol=1e-6` |
   | NMDS | stress within `0.02`; Procrustes correlation with R `> 0.95` | as stated |
   | Permutation tests (PERMANOVA) | test statistic elementwise; p-value within `0.02` at >= 9,999 permutations | as stated |
   | Rarefaction | invariants only: row sums == depth, dropped samples identical, no count exceeds original | exact |
   | DA methods | sign agreement and rank correlation of effect sizes; exact match only where the R method is deterministic | per method |

5. Any looser tolerance is written in the test with a one-line comment giving the reason.

# Why
R and NumPy random generators differ, so stochastic outputs can never match
bit-for-bit; demanding it would force skipping those tests. Comparing
invariants keeps them honest.

# Enforced by
- `tests/<module>/test_*_golden.py` files, marker `golden` (Phase 1, task 1.12).
- Missing golden file = test error, not skip.
