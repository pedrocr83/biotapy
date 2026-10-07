---
type: Module
title: da (differential abundance)
description: Four differential abundance methods, native LinDA and ANCOM-BC2 and the R bridges ALDEx2 and MaAsLin 3 (rpy2, extra `r`), that return one result table schema, and a consensus table counting where the methods agree; da writes no slot and filters nothing.
resource: /src/biotapy/da/
paths: ["src/biotapy/da/**"]
tags: [da, differential-abundance, linda, ancombc2, aldex2, maaslin3, rpy2]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-07T12:54:50Z }
commit: a1b54be
status: stable
---

# Responsibility

Owns the `bt.da.*` verbs of slices 3B and 3C:
- `linda` and `ancombc2` (native) and `aldex2` and `maaslin3` (R-only methods
  run through rpy2, extra `r`, [r-bridge-before-ports](/decisions/r-bridge-before-ports.md))
  test every feature between two groups (`linda`, `ancombc2` and `maaslin3`
  also against a numeric column) and return the same seven-column table, so
  any method's output can be compared with another's;
- `consensus` takes such tables and counts, per feature, how many methods call
  it and whether they agree
  ([da-consensus-agreement](/decisions/da-consensus-agreement.md)).

`da` writes no slot: results are returned, never stored
([data-model-slots](/contracts/data-model-slots.md), "DA results"). It does
not filter features or samples (run `pp.filter_features` once before any
method), takes no formula strings, and plots nothing ([pl](/modules/pl.md)'s
`consensus` draws the table). The bridges do not install or bundle R: the
user supplies R and its packages ALDEx2 and maaslin3 (`_r.py:r_function`).

# Entry points
- `_linda.py:linda` - least squares of each feature's log2 centred log-ratios
  on the design, with the effect's bias removed by the mode of all effects.
  Port of `MicrobiomeStat::linda`. `_linda.py:_mode` reproduces
  `modeest::mlv` (Gaussian mean shift) and `_linda.py:_shorth` its start value.
- `_ancombc.py:ancombc2` - a thin wrapper over scikit-bio's `ancombc2`; the
  group's row is taken from its result and divided by ln 2.
- `_aldex2.py:aldex2` - `ALDEx2::aldex` with the t-test; `effect` is its
  `diff.btw`, `pvalue` its `we.ep`, `se` all NaN. Two groups, no covariates.
  `_aldex2.py:_conditions` labels the group `"0"`/`"1"` and checks sizes.
- `_maaslin3.py:maaslin3` - `maaslin3::maaslin3`, abundance model only; takes
  the group's rows from the R return value. `_maaslin3.py:_string_levels`
  turns categories into strings.
- `_r.py:r_seed` (the integer for R's `set.seed`), `_r.py:r_function` (loads
  the R package and the R closure, returns a Python callable) and
  `_r.py:call_r` (converts pandas arguments and the data frame returned).
- `_consensus.py:consensus` - one row per feature in any table (first-seen
  order); `n_tested`, `n_significant`, `direction`, `consensus`, `conflict`.
- `_design.py:model` - validates `group`, `covariates` and `reference` once
  for every method and returns the frame and the `contrast` text.
  `_design.py:design_matrix` builds the matrix (`scale=True` for LinDA).
  `_design.py:dense_counts` checks counts and densifies `X` once; every method,
  native or bridged, goes through it.
- `_schema.py:result` builds the table; `_schema.py:validate_result` checks a
  table a user passes to `consensus`.

# Invariants
- **One schema.** Columns `effect`, `se`, `pvalue`, `qvalue`, `direction`,
  `method`, `contrast`, indexed by `feature`, rows in `var_names` order
  (`_schema.py:COLUMNS`, `_schema.py:result`).
- `effect` is a log2 fold change, with one method-specific exception:
  MaAsLin 3's is its coefficient minus the median coefficient
  (`_maaslin3.py:_MAASLIN`, `subtract_median = TRUE`), because that is what
  its p-value tests. `qvalue` is Benjamini-Hochberg over the finite p-values
  only, recomputed by biotapy for every method: not ALDEx2's `we.eBH`, not
  MaAsLin 3's `qval_individual` (`_schema.py:result`). `direction` is
  `sign(effect)`, 0 where `effect` is NaN.
- **NaN means not tested**: `effect`, `pvalue` and `qvalue` are NaN together,
  with `direction` 0. `_schema.py:validate_result` requires `effect` to be NaN
  exactly where `pvalue` is, and `qvalue` likewise. `se` is not checked: it is
  NaN for every feature of `aldex2` (ALDEx2 has none, `_aldex2.py:aldex2`).
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
  `dense_counts` raises for fewer than two features, for an empty sample and
  for repeated `var_names` or `obs_names`, since results are matched to
  features by name (`_design.py:dense_counts`, `_design.py:_require_unique`).
- `model` returns categorical columns **unordered**, without unused levels
  (`_design.py:_column`); see the ordered-factor gotcha.
- Numeric columns are scaled to unit variance by LinDA
  (`_design.py:design_matrix`, `scale=True` from `_linda.py:linda`) and by
  MaAsLin 3 itself (`_maaslin3.py:maaslin3` Notes); ANCOM-BC2 leaves them
  as they are, and `aldex2` refuses a numeric group
  (`_aldex2.py:_conditions`).
- The contrast reads `"<other level> vs <reference>"`, or the column name for
  a numeric group (`_design.py:_group`).
- `da` imports only `_core` and third-party packages
  ([module-boundaries](/contracts/module-boundaries.md)); `X` must be raw
  counts through `_core.require_counts` (`_design.py:dense_counts`).
- All four methods are compared with R in golden tests at measured
  tolerances; the seeded bridges elementwise
  ([r-golden-parity](/contracts/r-golden-parity.md)).
- **R bridge rules** (`_r.py`):
  - rpy2 is imported only through `import_optional(..., extra="r")` inside
    `_r.py:r_function` and `_r.py:call_r`, never at module level
    ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md));
    a missing R package becomes an `ImportError` naming the install line.
  - One `set.seed(r_seed(seed))` per call, inside the R closure
    (`_aldex2.py:_ALDEX`, `_maaslin3.py:_MAASLIN`); the wrapper restores the
    caller's `.Random.seed`, or its absence, on exit
    (`_r.py:_CATCH_WARNINGS`).
  - The converter is local (`_r.py:call_r`, `.context()`), never rpy2's
    global activation, so the user's own rpy2 session keeps its rules.
  - R warnings are caught as conditions and re-emitted as `UserWarning`; R
    messages are left alone. An R error is raised as `RuntimeError`, with the
    warnings raised before it in its message (`_r.py:r_function`).
  - ALDEx2 gets labels `"0"`/`"1"` (`_aldex2.py:_conditions`); MaAsLin 3
    gets string-level factors and plain names `x0`, `x1`, ... in its formula
    (`_maaslin3.py:maaslin3`), so an obs column named `body site` needs no
    quoting.
  - A feature MaAsLin 3 could not fit (error reported, or non-finite
    p-value) is untested (`_maaslin3.py:maaslin3`, `untested`); a feature
    ALDEx2 returns no row for (no read in any sample) is untested through the
    `reindex` (`_aldex2.py:aldex2`).

# Dependencies
- [core](/modules/core.md): `as_csr`, `require_counts`, `as_generator`,
  `import_optional`.
- scikit-bio (`ancombc2`), SciPy (`t`, `false_discovery_control`), NumPy and
  pandas; all core dependencies.
- rpy2 (extra `r`, GPL-2.0-or-later) and, installed by the user in R, the
  packages ALDEx2 and maaslin3; only `aldex2` and `maaslin3` need them
  ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).
- [da-consensus-agreement](/decisions/da-consensus-agreement.md) for what
  "agree" means; [pl](/modules/pl.md) draws `consensus`'s table.

# Verification

```bash
uv run --group test pytest tests/da src/biotapy/da -q
uv run --group test pytest -m "golden or network" tests/da -q
uv run --group test --extra r pytest -m r tests/da -q
```

The docs build also runs `da.linda`, `da.ancombc2` and `da.consensus` in
`docs/tutorials/differential_abundance.md`; each method's model, R settings
and measured agreement with R are written up in `docs/methods/`.
The second needs the pooch cache of GlobalPatterns (`BIOTAPY_DATA_DIR`); at
the 64fe39d check it gave `4 passed`. The third needs R with ALDEx2 and
maaslin3 installed (the `r-bridge` CI job, or the golden image in
`tests/r/`); there it gave `24 passed`. Without R the bridges' tests run
through the `fake_rpy2` fixture in the default run. Whether a faster
per-module command exists: unknown.

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
  LinDA and in MaAsLin 3 (both scale numeric columns), per unit in ANCOM-BC2;
  ALDEx2 refuses a numeric group. Do not rank or average them together.
  `_design.py:design_matrix`, `_maaslin3.py:maaslin3` Notes,
  `pl._consensus.py:consensus` Notes.
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
- **Ordered categoricals are unordered before R sees them.** R codes an
  ordered factor by polynomial contrasts, so a two-level ordered group
  reports an effect scaled by 1/sqrt(2) with no error. `_design.py:_column`
  calls `as_unordered()` for every categorical, which is why the
  pandas order of `group` never matters. Found in slice 3C's fix round for
  `maaslin3`; the native methods use category codes and were never affected.
- **R's random state is restored.** `set.seed` inside the closure would
  otherwise reseed the user's embedded R session; `_r.py:_CATCH_WARNINGS`
  puts `.Random.seed` back (or removes it) on exit. biotapy never touches
  NumPy's global state either (rules.md R3.4).
- **The seed integer is coupled to NumPy's generator stream.**
  `_r.py:r_seed` draws it with `integers(2**31 - 1)` from
  `as_generator(seed)`, and the goldens hard-code the result
  (`set.seed(1165433077)` for `seed=20260927`,
  `tests/da/test_aldex2_golden.py`, `tests/da/test_maaslin3_golden.py`). A
  NumPy change to that stream, or to the bound, fails the goldens loudly:
  recompute the integer and regenerate, do not loosen the test.
- **ALDEx2's reference swap is not exactly antisymmetric.** Its `diff.btw` resamples each
  group's pooled draws in label order (the draws and `we.ep` do not follow the
  labels), so swapping `reference` with the same seed
  moves effects by up to 0.65 log2 on GlobalPatterns' genera, where 15 of 636
  do not change direction; p-values and calls are the same (benchmark
  measurement of Checkpoint C). `_aldex2.py:aldex2` Notes.
- **MaAsLin 3's p-value moves about 1.2e-3 on a reference swap.** `effect`
  negates exactly, but the median test simulates around the coefficients, not
  their negatives (`pvalue` up to about 1.2e-3, `qvalue` up to about 3.3e-3 on the same data,
  calls unchanged). `_maaslin3.py:maaslin3` Notes.
- **ALDEx2 breaks on a factor `conds` and sorts labels by locale**, so
  `_aldex2.py:_conditions` passes the strings `"0"`/`"1"`; `aldex.effect`
  needs two samples in each level, which the same function checks.
- **MaAsLin 3 stops on `evaluate_only` with `warn_prevalence = TRUE`**, so the
  call sets it to FALSE. Its raw `coef` is not what its p-value tests (the
  median is subtracted), hence `subtract_median = TRUE`, and its
  `qval_individual` pools every covariate's p-values with the group's, hence
  biotapy's own BH over the group's. `_maaslin3.py:_MAASLIN`,
  `_maaslin3.py:maaslin3` Notes.
- **rpy2 sends non-string categories as sorted strings with only a
  `UserWarning`**, so a bool category would arrive in a different order than
  the one `model` set; `_maaslin3.py:_string_levels` renames categories to
  strings first, keeping the order.
- **A namespace that R loads with `::` during an rpy2 call prints "stack
  imbalance".** `_r.py:r_function` runs `importr(package)` before the call.
- **rpy2 without R's link headers falls back silently to an ABI mode that
  cannot load R**; CI installs the headers and sets `RPY2_CFFI_MODE=API` so
  the build fails instead (`.github/workflows/test.yaml`, job `r-bridge`).
- **LinDA densifies `X` and is native; the bridges densify too.**
  `_design.py:dense_counts` calls `toarray()` once (rpy2 has no sparse
  converter, so the bridges use it as well); the docstring gives the memory cost (about 5x the dense
  table at peak; ANCOM-BC2 4.5x to 7.5x, depending on the table's shape). rules.md R6.2 allows it: a native
  method whose algorithm needs the full table may densify once (the user
  approved that wording after Checkpoint B, 2026-10-05).
- **The tutorial quotes the four-method counts.** The docs build has no R, so
  `docs/tutorials/differential_abundance.md` prints the counts the `r` test
  `tests/da/test_consensus.py::test_four_methods_on_the_exit_gate_data`
  measures; both read the constants `FOUR_METHOD_*` in that file, and
  `test_the_tutorial_quotes_the_four_method_counts` fails until the tutorial
  matches them. A change to any method's numbers on GlobalPatterns updates
  the constants and the tutorial together.
- Replicate rows within a group give `se = 0` in LinDA and p-values that are
  floating-point noise; R does the same.

# Rejected
- Porting LinDA's paper-adaptive zero rule instead of the pseudocount path
  MicrobiomeStat 1.4 runs: it would differ from the R golden by up to 1.48
  log2 whenever a depth covariate exists (slice 3B decision 1 in
  [phase-3-stats](/roadmap/phase-3-stats.md)).
- Formula strings for the model: columns named by `group` and `covariates`
  keep validation in `_design.py:model`.
