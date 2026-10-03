---
type: Module
title: tl (tools)
description: Diversity, ordination and PERMANOVA over AnnData/TreeData - alpha, beta, UniFrac, PCoA, NMDS and PERMANOVA through scikit-bio and scikit-learn, returning results or writing the data-model-slots keys.
resource: /src/biotapy/tl/
paths: ["src/biotapy/tl/**"]
tags: [tl, diversity, ordination]
generated: { by: claude-code/claude-sonnet-5, at: 2026-10-03T13:00:00Z }
commit: 8f26269
status: stable
---

# Responsibility

Owns the `bt.tl.*` verbs that compute from a table or from a stored distance
matrix: alpha diversity, beta diversity and UniFrac (which write `obsp`),
ordination of an `obsp` matrix (PCoA, NMDS) and PERMANOVA over it. Owns no
reader and no transform: files are `io`'s, and anything that returns a
changed table (filter, rarefy, relative, tax_glom) is `pp`'s
([pp](/modules/pp.md)).

# Entry points

- `_alpha.py:alpha` - one value per sample and metric (`observed_features`,
  `shannon`, `simpson`, `chao1`, `faith_pd`) through
  `skbio.diversity.alpha_diversity`.
- `_beta.py:beta` - Bray-Curtis or binary Jaccard distances through
  `skbio.diversity.beta_diversity`.
- `_beta.py:unifrac` - unweighted or weighted UniFrac through the same call,
  with the tree from `_core.get_skbio_tree`.
- `_ordination.py:pcoa` - `skbio.stats.ordination.pcoa` of `obsp[distance]`.
- `_ordination.py:nmds` - non-metric SMACOF (`sklearn.manifold.MDS`) of
  `obsp[distance]`.
- `_permanova.py:permanova` - `skbio.stats.distance.permanova` of
  `obsp[distance]` against an `obs` column.
- `_beta.py:stored_distances` - private, shared by `pcoa`, `nmds` and
  `permanova` (same subpackage, [module-boundaries](/contracts/module-boundaries.md)):
  reads `obsp[key]` as a `DistanceMatrix`, raising `KeyError` that names the
  `bt.tl...` call writing a missing key (`_beta.py:_WRITTEN_BY`) and
  `ValueError` on NaN distances.

# Invariants

- Everything defaults to `inplace=False` and returns its result
  ([pure-by-default](/decisions/pure-by-default.md)). With `inplace=True` a
  function writes only the [data-model-slots](/contracts/data-model-slots.md)
  keys and returns `None`: `obs["alpha_<metric>"]`; `obsp["braycurtis" |
  "jaccard" | "unweighted_unifrac" | "weighted_unifrac"]` (`_beta.py:_store`);
  `obsm["X_pcoa"]` plus `uns["biotapy"]["pcoa"]`; `obsm["X_nmds"]` plus
  `uns["biotapy"]["nmds"]`. The `assert_unchanged` fixture
  (`tests/conftest.py`) checks inputs, including the key set of
  `uns["biotapy"]`.
- scikit-bio only ever gets dense input (rules.md R6.2): `alpha` densifies
  rows in chunks of at most `2**20` values (`_alpha.py:_CHUNK_VALUES`);
  `beta` and `unifrac` densify `X` once, since pairwise distances need every
  row.
- Trees come only from `_core.get_skbio_tree`; nothing in `tl` imports
  `treedata` or `networkx` ([tree-access](/contracts/tree-access.md)).
- `pcoa` asks scikit-bio for at most `n_obs - 1` axes, so
  `proportion_explained` divides by the trace of the centred matrix, like
  `ape::pcoa`'s `Relative_eig`; with all `n` axes scikit-bio would divide by
  the sum of the positive eigenvalues instead. `_ordination.py:pcoa`
- `nmds` keeps the best of 20 random starts (`_ordination.py:_NMDS_STARTS`,
  vegan's `try = 20`) and passes every `MDS` argument explicitly, because
  scikit-learn 1.9/1.10 change the defaults of `n_init` and `init`.
- `permanova` has no `inplace`: it returns a test result, not per-sample or
  per-pair values.
- `faith_pd` runs on presence/absence (`(dense > 0)` as int64, per chunk),
  which is exact because Faith PD depends on presence only, so it runs on
  any abundance. `_alpha.py:alpha`
- Weighted UniFrac and the count metrics (`observed_features`, `chao1`) need
  whole-number counts: they call `_core.require_counts`, which checks both
  `x_kind` and the values of `X`. Unweighted UniFrac and Jaccard are
  qualified to presence by scikit-bio and need no counts.

# Dependencies

- [core](/modules/core.md): `as_csr`, `require_counts`, `get_skbio_tree`,
  `as_generator`, `TreeData`.
- scikit-bio: `alpha_diversity`, `beta_diversity`, `DistanceMatrix`, `pcoa`,
  `permanova`.
- scikit-learn: `sklearn.manifold.MDS`, with an `int` seed drawn from
  `as_generator`, because its `random_state` rejects a `np.random.Generator`.

# Verification

`uv run --group test pytest tests/tl -q`, and the R golden tests with a pooch
cache: `BIOTAPY_DATA_DIR=<cache> uv run --group test pytest -m golden tests/tl -q`.

# Gotchas

- All-zero samples give scikit-bio's values, with no custom mapping: two of
  them are NaN apart under Bray-Curtis and 0 apart under Jaccard and both
  UniFracs; `alpha` gives 0 for `observed_features`, `chao1` and `faith_pd`
  and NaN for `shannon` and `simpson` (phyloseq: Shannon 0, Simpson 1).
  NaN distances then make `pcoa`, `nmds` and `permanova` raise
  (`_beta.py:stored_distances`).
- PCoA axis signs are arbitrary, in R too; tests align them per axis before
  comparing.
- scikit-bio sets negative PCoA eigenvalues, and their coordinates, to 0,
  whereas `ape::pcoa` reports those eigenvalues as negative.
- NMDS on `toy()` is degenerate (stress about 0.002), so unit tests check
  NMDS on planar points (`tests/tl/test_ordination.py:test_nmds_recovers_a_planar_configuration`).
- The PERMANOVA golden test takes 10-20 s at 9,999 permutations (18 s at
  Checkpoint C).
- `faith_pd` on 5,000 x 50,000 takes about 48 s (asv, Task 1.21). scikit-bio
  re-indexes and re-validates the tree on every `alpha_diversity` call, once per
  chunk of `2**20 // n_vars` samples, and that is about 97% of the time;
  converting the tree once (`get_skbio_tree`) takes 0.46 s. A fix needs a
  profile-driven perf task (rules.md R10.1), not a change here.
- scikit-bio 0.7.4's tree code, `skbio.diversity._phylogenetic._nodes_by_counts`,
  casts abundances to int64. Faith PD and both UniFrac engines (Cython and
  numba) share it, so fractions are truncated silently. This is why
  `faith_pd` gets presence/absence and weighted UniFrac calls
  `require_counts` (Checkpoint C fix C1).
- `obs["alpha_*"]` columns survive `pp.filter_features`, `pp.rarefy` and
  `pp.tax_glom`, as mia keeps `colData`; `_core.feature_subset` keeps `obs`.
  They then describe the old features and must be recomputed.
- A numeric `obs[grouping]` (not bool) makes `permanova` raise `TypeError`:
  it would become one group per value, while `adonis2` fits a numeric term as
  continuous. `.astype("category")` gives groups. `_permanova.py:permanova`
