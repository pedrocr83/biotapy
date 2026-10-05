---
type: Contract
title: R golden parity
description: Every computation with an R equivalent is tested against gzip CSV golden files exported from pinned R; deterministic outputs match numerically, stochastic outputs match invariants; HUMAnN-parity functions are tested the same way against files exported from pinned HUMAnN 3.9.
tags: [testing, r, validation]
status: stable
paths: ["tests/r/**", "tests/golden/**", "tests/**/test_*.py"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T15:08:33Z }
commit: 2216894
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
   `Biostrings`, `vegan` and `ape`, plus CRAN `picante` for Faith PD and
   Bioconductor `philr` for `pp.philr`, with the `libuv1` runtime library its
   `fs` binary loads); a new golden function that needs another package adds
   it in its own commit (rules.md R2.3).
1b. HUMAnN golden files (`fn.func_glom`, `fn.renorm`) are produced by
   `tests/humann/export_golden.py`, run with
   `uv run --no-project --with humann==3.9 --with pandas==3.0.6`, never
   in CI. HUMAnN's utility scripts are pure Python and need no
   database, so there is no container. Output:
   `tests/golden/humann/<name>.csv.gz` (samples as rows, row ids without
   their `": name"`) and `tests/golden/humann/VERSIONS.txt`.
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
   | Deterministic numeric (glom sums, relative, CLR, alpha, Bray-Curtis, UniFrac) | elementwise | `rtol=1e-7` |
   | PCoA coordinates | per axis, up to sign flip; eigenvalues elementwise | `rtol=1e-6` |
   | PhILR balances | matched by partition (the taxa in each numerator and denominator, so signs must agree too), then elementwise | `rtol=1e-7`; `atol=1e-12`, because a balance between absent taxa is 0 here and about 1e-16 in R |
   | NMDS | stress within `0.02`; Procrustes correlation with R `> 0.95` | as stated |
   | Permutation tests (PERMANOVA) | test statistic elementwise; p-value within `0.02` at >= 9,999 permutations | as stated |
   | Rarefaction | invariants only: row sums == depth, dropped samples identical, no count exceeds original | exact |
   | HUMAnN parity (func_glom, renorm) | elementwise, matched by row id | rtol=1e-7; renorm rtol=5e-6, because humann_renorm_table prints %.6g |
   | DA methods | sign agreement and rank correlation of effect sizes; exact match only where the R method is deterministic | per method |

5. Any looser tolerance is written in the test with a one-line comment giving the reason.
6. Golden files hold numbers derived from third-party example data, never the
   raw files: those are fetched at test time by pooch. Derived numbers may be
   complete for a dataset (`global_patterns/relative.csv.gz` holds every
   nonzero proportion of GlobalPatterns, AGPL-3 via phyloseq); the user
   accepted this for this BSD-3 repository on 2026-09-27. Test fixtures under `tests/data/`
   stay synthetic, except small files copied under a permissive licence
   with a `NOTICE.txt` beside them (`tests/data/humann`: HUMAnN's MIT test
   data; `tests/data/metaphlan`: a MetaPhlAn 4.0.6 profile from HUMAnN's MIT
   test data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0). PICRUSt2
   (GPL-3) fixtures are always synthetic, written from its documented
   column headers.
7. `pl` functions have an R equivalent but no golden test. They draw numbers
   that `tl` stores, and `tl`'s golden tests check those numbers (controller
   ruling 2026-09-27; rules.md R11.2). `pl.contributions` has no R
   equivalent and draws `fn.contributions`, which has none either (hand-computed
   cases and Hypothesis properties check it, `tests/fn/test_contributions.py`).
8. A reader's R parity may come from a tool's own golden files downstream,
   when an R golden for the reader would only re-check parsed numbers.
   `io.read_humann` (R equivalent `mia::importHUMAnN`) is checked this way:
   the HUMAnN golden tests of `fn.func_glom` and `fn.renorm` read their
   inputs through it and compare with HUMAnN's own output. mia is not added
   to the R image (user-approved 2026-10-03). `io.read_metaphlan` (R equivalent
   `mia::importMetaPhlAn`) is checked against MetaPhlAn's own output: its
   leaves, grouped by `pp.tax_glom` to each rank, equal the clade rows the
   profile prints (`tests/io/test_metaphlan.py`, atol `1e-6` because
   MetaPhlAn rounds each percentage to 5 decimals). `io.read_picrust2` and
   `io.read_picrust2_traits` have no R equivalent; invariants on synthetic
   files check them.

# Why
R and NumPy random generators differ, so stochastic outputs can never match
bit-for-bit; demanding it would force skipping those tests. Comparing
invariants keeps them honest.

# Enforced by
- `tests/<module>/test_*_golden.py` files, marker `golden` (Phase 1, task 1.12).
  They also carry `network`, because the golden inputs (e.g. GlobalPatterns)
  are downloaded by pooch; they run in the network CI job
  (`pytest -m "network or golden"`).
- `tests/fn/*_golden.py`, marker `golden` only: their inputs are committed, so
  they run in every CI job.
- Missing golden file = test error, not skip.
