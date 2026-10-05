---
type: Module
title: da (differential abundance)
description: Two native differential abundance methods, LinDA and ANCOM-BC2, that return one result table schema, and a consensus table counting where the methods agree; da writes no slot and filters nothing.
resource: /src/biotapy/da/
paths: ["src/biotapy/da/**"]
tags: [da, differential-abundance, linda, ancombc2]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-05T22:20:19Z }
commit: 711b643
status: stable
---

# Responsibility

Owns the `bt.da.*` verbs of slice 3B:
- `linda` and `ancombc2` test every feature between two groups (or against a
  numeric column) and return the same seven-column table, so any method's
  output can be compared with another's;
- `consensus` takes such tables and counts, per feature, how many methods call
  it and whether they agree
  ([da-consensus-agreement](/decisions/da-consensus-agreement.md)).

`da` writes no slot: results are returned, never stored
([data-model-slots](/contracts/data-model-slots.md), "DA results"). It does
not filter features or samples (run `pp.filter_features` once before any
method), takes no formula strings, and plots nothing ([pl](/modules/pl.md)'s
`consensus` draws the table). The R bridges (`aldex2`, `maaslin3`) are later
tasks; this concept is updated for them in task 3.14.

# Entry points
- `_linda.py:linda` - least squares of each feature's log2 centred log-ratios
  on the design, with the effect's bias removed by the mode of all effects.
  Port of `MicrobiomeStat::linda`. `_linda.py:_mode` reproduces
  `modeest::mlv` (Gaussian mean shift) and `_linda.py:_shorth` its start value.
- `_ancombc.py:ancombc2` - a thin wrapper over scikit-bio's `ancombc2`; the
  group's row is taken from its result and divided by ln 2.
- `_consensus.py:consensus` - one row per feature in any table (first-seen
  order); `n_tested`, `n_significant`, `direction`, `consensus`, `conflict`.
- `_design.py:model` - validates `group`, `covariates` and `reference` once
  for every method and returns the frame and the `contrast` text.
  `_design.py:design_matrix` builds the matrix (`scale=True` for LinDA).
  `_design.py:dense_counts` checks counts and densifies `X` once.
- `_schema.py:result` builds the table; `_schema.py:validate_result` checks a
  table a user passes to `consensus`.

# Invariants
- **One schema.** Columns `effect`, `se`, `pvalue`, `qvalue`, `direction`,
  `method`, `contrast`, indexed by `feature`, rows in `var_names` order
  (`_schema.py:COLUMNS`, `_schema.py:result`).
- `effect` is a log2 fold change. `qvalue` is Benjamini-Hochberg over the
  finite p-values only (`_schema.py:result`). `direction` is `sign(effect)`,
  0 where `effect` is NaN.
- **NaN means not tested**, in all four float columns together with
  `direction` 0. `_schema.py:validate_result` requires `effect` to be NaN
  exactly where `pvalue` is, and `qvalue` likewise.
- No feature is dropped: an untested feature keeps its row.
  `tests/da/test_schema.py`.
- Calls are strict: `qvalue < alpha`, never `<=`
  (`_consensus.py:consensus`).
- A group is categorical with exactly two levels, or numeric. Missing values
  in a used `obs` column, a constant column, a repeated column, fewer samples
  than model terms and collinear columns all raise
  (`_design.py:model`).
- `model` raises `TypeError` for `covariates` that is not a list of column
  names and for a non-string `reference`; `ValueError` for a constant column
  and for a `reference` that is not a level or is given for a numeric group
  (`_design.py:model`, `_design.py:_check_varies`, `_design.py:_group`).
  `dense_counts` raises for fewer than two features
  (`_design.py:dense_counts`).
- Numeric columns are scaled to unit variance in LinDA only
  (`_design.py:design_matrix`, `scale=True` from `_linda.py:linda`).
- The contrast reads `"<other level> vs <reference>"`, or the column name for
  a numeric group (`_design.py:_group`).
- `da` imports only `_core` and third-party packages
  ([module-boundaries](/contracts/module-boundaries.md)); `X` must be raw
  counts through `_core.require_counts` (`_design.py:dense_counts`).
- Both methods are compared with R in golden tests at measured tolerances
  ([r-golden-parity](/contracts/r-golden-parity.md)).

# Dependencies
- [core](/modules/core.md): `as_csr`, `require_counts`.
- scikit-bio (`ancombc2`), SciPy (`t`, `false_discovery_control`), NumPy and
  pandas; all core dependencies.
- [da-consensus-agreement](/decisions/da-consensus-agreement.md) for what
  "agree" means; [pl](/modules/pl.md) draws `consensus`'s table.

# Verification

```bash
uv run --group test pytest tests/da src/biotapy/da -q
uv run --group test pytest -m "golden or network" tests/da -q
```

The second needs the pooch cache of GlobalPatterns (`BIOTAPY_DATA_DIR`);
at the 927e5ae check it gave `4 passed`. Whether a faster per-module command
exists: unknown.

# Gotchas
- **ANCOM-BC2 is not antisymmetric in `reference`.** Swapping it changes more
  than the sign: on GlobalPatterns' genera (`host`) the calls at q < 0.05 go
  from 208 to 230, and the effect and its swapped effect sum to -0.40 to -0.37
  log2 (0 if only the sign changed). R does the same, so it is inherited, not a
  bug; consensus inherits it (104 against 112 calls with LinDA).
  `_ancombc.py:ancombc2` Notes. Fix the reference on the biology first.
- **No residual degrees of freedom.** scikit-bio reports such a feature with a
  finite effect and NaN p; `ancombc2` makes it all-NaN (untested) so the
  schema's NaN rule holds. scikit-bio reports an unfitted feature with
  p = 1; that too becomes untested, left out of the BH correction, where R
  counts it. `_ancombc.py:ancombc2` (`untested`).
- **A numeric group's effect has different units.** Per standard deviation in
  LinDA (it scales numeric columns), per unit in ANCOM-BC2. Do not rank or
  average them together. `_design.py:design_matrix`, `pl._consensus.py:consensus`
  Notes.
- **scikit-bio fit errors are re-raised** as `ValueError` with the prefix
  `da.ancombc2:` and `from err`. `_ancombc.py:ancombc2`.
- **ANCOM-BC2's E-M is capped at 100 iterations**, R's default, and has not
  converged on some data; with 1,000 iterations R moves every effect by about
  -0.28 log2 on GlobalPatterns. Kept for parity. `_ancombc.py:ancombc2` Notes.
- **MicrobiomeStat's adaptive zero imputation is unreachable** (it compares
  `"Imputation"` with `"imputation"`), so 0.5 is added whenever `X` holds a
  zero; biotapy computes what R returns. `_linda.py:linda`.
- **`modeest` quirks are ported**: the loop returns the previous mean-shift
  value, and tied shorth starts are averaged then truncated.
  `_linda.py:_mode`, `_linda.py:_shorth`.
- **patsy puts numeric terms after categorical ones** in scikit-bio's design
  (noted in the slice 3B design, [phase-3-stats](/roadmap/phase-3-stats.md)).
  `ancombc2` renames columns to `x0`, `x1`, ... and picks the group's rows by
  name (`x0` or `x0[...]`), never by position. `_ancombc.py:ancombc2`.
- **LinDA densifies `X` and is native.** `_design.py:dense_counts` calls
  `toarray()` once; the docstring gives the memory cost (about 5x the dense
  table at peak; ANCOM-BC2 about 7x). rules.md R6.2 allows it: a native
  method whose algorithm needs the full table may densify once (the user
  approved that wording after Checkpoint B, 2026-10-05).
- Replicate rows within a group give `se = 0` in LinDA and p-values that are
  floating-point noise; R does the same.

# Rejected
- Porting LinDA's paper-adaptive zero rule instead of the pseudocount path
  MicrobiomeStat 1.4 runs: it would differ from the R golden by up to 1.48
  log2 whenever a depth covariate exists (slice 3B decision 1 in
  [phase-3-stats](/roadmap/phase-3-stats.md)).
- Formula strings for the model: columns named by `group` and `covariates`
  keep validation in `_design.py:model`.
