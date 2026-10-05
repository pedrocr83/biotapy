---
type: Module
title: pp (preprocessing)
description: Pure transforms over AnnData/TreeData that scale abundances per sample, take the centred log-ratio or the phylogenetic balances (PhILR), aggregate features along the taxonomy, or filter samples/features.
resource: /src/biotapy/pp/
paths: ["src/biotapy/pp/**"]
tags: [pp]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T15:43:25Z }
commit: 6ade269
status: stable
---

# Responsibility

Owns the `bt.pp.*` verbs that transform an AnnData/TreeData's abundance
table: `relative` (adds a layer, keeps every feature), `clr` (adds a dense
layer), `philr` (adds an `obsm` embedding over the tree), `tax_glom`
(aggregates features to a taxonomic rank, drops derived slots). Owns
filtering (`filter_features`, `filter_samples`) and rarefaction (`rarefy`);
does NOT own any diversity/ordination computation ([tl](/modules/tl.md)).

# Entry points

- `_transform.py:relative` - per-sample scaling to relative abundance (each
  stored value divided by its float64 sample total through
  `_core._matrix.py:divide_rows`, never multiplied by a reciprocal, which
  overflows for a subnormal total); all-zero samples stay all-zero, which differs from phyloseq's
  `transform_sample_counts` (returns `NaN` there).
- `_transform.py:clr` - centred log-ratio into `layers["clr"]`, through
  scikit-bio's `clr` on `_transform.py:pseudocounted`.
- `_transform.py:pseudocounted` - the check-and-densify step `clr` and `philr`
  share: validates `pseudocount`, rejects negative or non-finite `X`, adds the
  pseudocount to every value and returns one dense float64 array; `columns=`
  reorders features while still sparse.
- `_philr.py:philr` - PhILR balances into `obsm["X_philr"]`, through
  scikit-bio's `tree_basis` on the tree from `_core.get_skbio_tree`.
- `_philr.py:_binary_tree` - copy of the tree without one-child nodes, children
  in order; raises on a node with more than two.
- `_glom.py:tax_glom` - aggregate features sharing a lineage down to `rank`
  into their most abundant member.
- `_filter.py:filter_features` - prevalence and total thresholds, through
  `feature_subset`.
- `_filter.py:filter_samples` - a depth threshold, through AnnData indexing,
  keeping every slot.
- `_rarefy.py:rarefy` - subsample every sample to the same depth, through
  `skbio.stats.subsample_counts`, then `feature_subset`.

# Invariants

- Every function is pure ([pure-by-default](/decisions/pure-by-default.md)):
  returns a new object, input unchanged. Checked via the `assert_unchanged`
  fixture (`tests/conftest.py`) in `tests/pp/test_{transform,glom,filter,rarefy}.py`.
- `tax_glom` groups by the full lineage string joined down to `rank`, not by
  the rank label alone - phyloseq semantics: `";_;".join`, missing value ->
  `"NA"` - so same-named taxa in different lineages (e.g. two unrelated
  "uncultured" genera) stay separate. `_glom.py:tax_glom`
- The group archetype is the most abundant member by column sum across all
  samples, first on ties; kept features keep their original column order.
  `_glom.py:tax_glom`, pinned against a naive reference in
  `tests/pp/test_glom.py:test_matches_naive_phyloseq_reference` and
  `tests/pp/test_glom.py:test_interleaved_taxa_columns_end_up_in_original_order`.
- `tax_glom` goes through `_core.feature_subset`, so a TreeData's tree is
  pruned to the archetypes' leaves plus ancestors rather than rebuilt - see
  [core](/modules/core.md) for what that keeps.
- `filter_features`/`filter_samples` thresholds are inclusive.
  `_filter.py:filter_features`, `_filter.py:filter_samples`
- `filter_features` counts prevalence from stored non-zero values and divides
  by `n_obs`, rather than multiplying the threshold by `n_obs`, because a
  fraction like `3 / 10` is not exactly representable and multiplying can put
  the boundary on the wrong side. `_filter.py:filter_features`,
  `tests/pp/test_filter.py:test_decimal_prevalence_boundary_is_kept`.
- Nothing passing either filter raises `ValueError`; neither ever returns an
  empty object. `_filter.py:filter_features`, `_filter.py:filter_samples`
- `rarefy` subsamples without replacement only (phyloseq defaults to
  `replace = TRUE`); `depth` defaults to the smallest non-zero sample depth,
  so all-zero samples are skipped rather than rarefying everything to 0.
  Samples with strictly fewer than `depth` reads are dropped with one
  `UserWarning` naming up to 5 of them; a sample at exactly `depth` is kept.
  Features left all-zero after subsampling are dropped through
  `feature_subset`, and `X` becomes `int64` counts. `_rarefy.py:rarefy`

- `clr` and `philr` add the pseudocount to every value of `X`, zeros or not
  (equal to `vegan::decostand` and `philr::philr`, pinned by the golden tests
  `tests/pp/test_transform_golden.py` and `tests/pp/test_philr_golden.py`). A
  pseudocount of 0 raises when `X` holds a zero. A pseudocount above the
  smallest non-zero value of `X` gives one `UserWarning` through `warn_user`
  (the default 0.5 meets relative abundances). `_transform.py:pseudocounted`
- `layers["clr"]` is dense float64, because CLR has no zeros; an all-zero
  sample gives an all-zero row. `_transform.py:clr`
- `obsm["X_philr"]` is a samples x balances `DataFrame`, one column per
  internal node with two children, in preorder, named by the node names of
  `vart["phylo"]`. A balance is positive when the node's first child is the
  more abundant; this is `philr::philr`'s sign. `_philr.py:philr`
- `philr` skips one-child nodes (a subsetted TreeData keeps them; R's
  `drop.tip` removes them) and raises `ValueError` for a node with more than
  two children, the root included, naming up to 3 and the total. It does not
  resolve polytomies. `_philr.py:_binary_tree`,
  `_core/_tree.py:get_skbio_tree` (called with `split_root=False`)

# Dependencies

- [core](/modules/core.md): `as_csr`, `sum_by`, `argmax_by`, `split_ranks`,
  `feature_subset`, `add_provenance`, `as_generator`, `require_counts`,
  `warn_user`, `divide_rows`, `get_skbio_tree`.

# Verification

`uv run --group test pytest tests/pp -q` plus `uvx prek run --all-files`.

# Gotchas

- `relative`'s all-zero-sample behavior is a deliberate departure from
  phyloseq, stated in the docstring `Notes` - do not "fix" it to return NaN.
- `tax_glom`'s output column order follows the *archetypes'* original
  position in `var_names`, not the order groups are first encountered while
  scanning lineages; `bt.datasets.toy()`'s taxa are contiguous per phylum, so
  `tests/pp/test_glom.py` builds its own interleaved fixture to cover this.
- `filter_samples` keeps every slot, so a kept `obsm` ordination and its
  `uns["biotapy"]["pcoa"|"nmds"]` summary were computed with the dropped
  samples included; recompute them (`_filter.py:filter_samples` docstring).
- **`tree_basis` and `TreeNode.prune` fight PhILR's conventions.** scikit-bio's
  `tree_basis` puts a node's first child in the denominator, `philr::philr` in
  the numerator, so `philr` negates the result. `TreeNode.prune` would drop
  one-child nodes but moves the surviving child to the end of its parent's
  list, which would break the first-child sign rule, so `_binary_tree`
  rebuilds the tree instead.
  `_philr.py:philr`, `_philr.py:_binary_tree`.
- `bt.datasets.toy()`'s root has three children, so `pp.philr` raises on it;
  the docstring example slices to `f1`-`f6`, which sit under two.
  `_philr.py:philr`.
- **Plain slicing keeps stale derived slots.** `tdata[:, cols]` keeps
  `layers["clr"]` and `obsm["X_philr"]`, which no longer match the features;
  only `filter_features` and the other `feature_subset` callers drop them.
  Rerun `clr`/`philr` after slicing. `_core/_slots.py:feature_subset`,
  `tests/pp/test_philr.py` (the `filter_features` case).
- PhILR balance names are `vart["phylo"]`'s node names, not R's
  `makeNodeLabel` names (`n1`, `n2`, ...); the golden test matches balances
  by the partition of tips they split, not by name
  ([r-golden-parity](/contracts/r-golden-parity.md)). To binarise a tree in
  Python use scikit-bio's `TreeNode.bifurcate()`; R users use
  `ape::multi2di`. `_philr.py:philr`.
- **Peak memory is several dense arrays.** `clr` peaks at 3.1x to 4.8x one
  dense `X` (measured on 400 x 512 and 2,000 x 2,000; the most when `X` is
  mostly non-zero). `philr` peaks at about five arrays plus scikit-bio's
  sparse basis, which holds one value per tip under each node (113 MB on
  GlobalPatterns' 26 x 19,216). Both state it in their `Notes`.
  `_transform.py:clr`, `_philr.py:philr`.
- `pseudocounted` lives in `_transform.py` and `_philr.py` imports it from
  there, which module-boundaries allows inside one subpackage. If slice 3B's
  `da.linda` reuses it, it moves to `_core` (two or more subpackages,
  [module-boundaries](/contracts/module-boundaries.md) rule 2).
- `filter_features` keeps a feature when `present / n_obs >= min_prevalence`;
  phyloseq's `sum(x > 0) >= p * length(x)` can drop it at an exact boundary
  through floating point (7 of 25 samples at `p = 0.28`).
