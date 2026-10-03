---
type: Module
title: pl (plots)
description: Plots of what tl and pp stored - stacked bars, heatmap, richness, ordination and scree - drawn with matplotlib on the given or a new Axes, computing nothing.
resource: /src/biotapy/pl/
paths: ["src/biotapy/pl/**"]
tags: [pl, plots, matplotlib]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-03T08:10:00Z }
commit: 2b9fc24
status: stable
---

# Responsibility

Owns the `bt.pl.*` verbs that draw stored results:
- `bar` and `heatmap` draw `X` or a layer;
- `richness` draws `obs["alpha_<metric>"]`;
- `ordination` draws `obsm["X_pcoa" | "X_nmds"]` with its `uns["biotapy"]`
  summary;
- `scree` draws `uns["biotapy"]["pcoa"]`.

It computes no diversity, distance or ordination; those are
[tl](/modules/tl.md)'s. A missing slot is a `KeyError` naming the `tl` or
`pp` call that writes it; for a `pp` layer it names the assignment
(`adata = bt.pp.relative(adata)`), because `pp` returns a new object.

# Entry points
- `_abundance.py:bar` - stacked bars. `_segments` sums features per `fill`
  group (a `var` or an `obs` column); `_by_x` sums samples per `x` group.
- `_abundance.py:heatmap` - `imshow` of the table densified once, features x
  samples, on a `LogNorm` scale with phyloseq's colours.
- `_richness.py:richness` - one point per sample; NaN values are left out,
  and a `color` group left with no point gets no legend entry
  (`_common.py:scatter` skips empty groups).
- `_ordination.py:ordination` and `_ordination.py:scree` - both read through
  `_ordination.py:_stored`.
- `_common.py` - helpers the topic files share: `new_axes`, `table`,
  `obs_groups`/`groups`, `scatter` and `label_ticks`.

# Invariants
- `import biotapy` never imports matplotlib. Annotations use a
  `TYPE_CHECKING` import, and functions import matplotlib when they draw
  (`tests/pl/test_init.py`).
- Every function takes `ax=None` last and keyword-only, and returns the
  `Axes`. With `ax=None` it draws on a new pyplot figure
  (`_common.py:new_axes`), which notebooks display.
- Nothing is written to the AnnData
  ([pure-by-default](/decisions/pure-by-default.md)); each function has a test
  asserting the input unchanged (`assert_unchanged`, `tests/conftest.py`).
- Groups (`_common.py:groups`) keep a categorical's order and sort other
  values. Missing values are a last `NA` group, drawn grey. A numeric column
  raises `TypeError` through `_core.require_categorical`, as `tl.permanova`
  does.
- `bar` heights equal the sample (or `x` group) totals of the plotted table:
  features with a missing rank are a group, not dropped (a Hypothesis test).
- `heatmap` keeps `obs`/`var` order and densifies the table once (rules.md
  R6.2); it raises `ValueError` naming `adata` (for `X`) or `layer=` when
  nothing is positive.
- Past 250 names an axis gets no tick labels (`_common.py:MAX_LABELS`,
  phyloseq's `max.label`).

# Dependencies
- [core](/modules/core.md): `as_csr`, `sum_by`, `require_categorical`.
- matplotlib `>=3.8`, a runtime dependency
  ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).
- The slots [tl](/modules/tl.md) and `pp.relative` write
  ([data-model-slots](/contracts/data-model-slots.md)).

# Verification

`uv run --group test pytest tests/pl src/biotapy/pl -q`. The tutorials draw
every function on GlobalPatterns, enterotype and esophagus in the docs build:
`BIOTAPY_DATA_DIR=<cache> uv run --group doc sphinx-build -W -b html docs docs/_build/html`.

# Gotchas
- `ax.bar` adds one `Rectangle` per bar, one group after another, so tests
  reshape the `ax.patches` heights to bars x groups.
- `LogNorm` masks zeros, which the colormap draws in its "bad" colour
  (black). On an all-zero table matplotlib would fail inside `colorbar` with
  `Invalid vmin or vmax`, so `heatmap` raises its own `ValueError` first.
- `ax.boxplot` on the same axes replaces the group tick labels unless it
  gets `manage_ticks=False` (the vignette's richness cell).
- `bar` with hundreds of `fill` groups is slow: GlobalPatterns by genus is
  984 groups and took 15.6 s on a loaded machine (about 8 s earlier). Aggregate first.
- matplotlib 3.8.0-3.8.3 pin `numpy<2`, so with NumPy 2 the resolver takes
  3.8.4 or later.
- `pl` has R equivalents but no golden tests: it draws numbers that `tl`'s
  golden tests check ([r-golden-parity](/contracts/r-golden-parity.md)
  statement 7).
