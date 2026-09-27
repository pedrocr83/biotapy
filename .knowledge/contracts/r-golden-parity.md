---
type: Contract
title: R golden parity
description: Every function with an R equivalent is tested against gzip CSV golden files exported from pinned R; deterministic outputs match numerically, stochastic outputs match invariants.
tags: [testing, r, validation]
status: stable
paths: ["tests/r/**", "tests/golden/**", "tests/**/test_*.py"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-27T09:42:00Z }
commit: 625429c
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Statement
1. Golden files are produced by `tests/r/export_golden.R`, run only inside the
   pinned container `tests/r/Dockerfile`, never in normal CI. The image
   installs only what the current golden files need (`phyloseq`, which brings
   `Biostrings`); a new golden function that needs another package adds it in
   its own commit (rules.md R2.3).
2. Output: `tests/golden/<dataset>/<function>.csv.gz`, samples as rows
   (see [samples-as-rows](/decisions/samples-as-rows.md)), plus
   `tests/golden/VERSIONS.txt` listing R and package versions. Golden files
   are gzip CSV, not parquet: the user declined `pyarrow` (2026-09-27), so the
   R image needs no `arrow` and tests read the files with plain pandas.
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
  They also carry `network`, because the golden inputs (e.g. GlobalPatterns)
  are downloaded by pooch; they run in the network CI job
  (`pytest -m "network or golden"`).
- Missing golden file = test error, not skip.
