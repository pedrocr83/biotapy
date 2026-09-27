---
type: Module
title: pp (preprocessing)
description: Pure transforms over AnnData/TreeData that scale abundances per sample, aggregate features along the taxonomy, or filter samples/features.
resource: /src/biotapy/pp/
paths: ["src/biotapy/pp/**"]
tags: [pp]
generated: { by: claude-code/claude-sonnet-5, at: 2026-09-27T17:46:59Z }
commit: d28af22
status: stable
---

# Responsibility

Owns the `bt.pp.*` verbs that transform an AnnData/TreeData's abundance
table: `relative` (adds a layer, keeps every feature), `tax_glom`
(aggregates features to a taxonomic rank, drops derived slots). Owns
filtering (`filter_features`, `filter_samples`) and rarefaction (`rarefy`);
does NOT own any diversity/ordination computation (`tl`, Slice 1C).

# Entry points

- `_transform.py:relative` - per-sample scaling to relative abundance;
  all-zero samples stay all-zero, which differs from phyloseq's
  `transform_sample_counts` (returns `NaN` there).
- `_glom.py:tax_glom` - aggregate features sharing a lineage down to `rank`
  into their most abundant member.
- `_filter.py:filter_features` - prevalence and total thresholds, through
  `feature_subset`.
- `_filter.py:filter_samples` - a depth threshold, through AnnData indexing,
  keeping every slot.
- `_rarefy.py:rarefy` - subsample every sample to the same depth, through
  `skbio.stats.subsample_counts`, then `feature_subset`.

# Invariants

- Both functions are pure ([pure-by-default](/decisions/pure-by-default.md)):
  return a new object, input unchanged. Checked via the `assert_unchanged`
  fixture (`tests/conftest.py`) in `tests/pp/test_transform.py` and
  `tests/pp/test_glom.py`.
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

# Dependencies

- [core](/modules/core.md): `as_csr`, `sum_by`, `argmax_by`, `split_ranks`,
  `feature_subset`, `add_provenance`, `as_generator`, `require_counts`,
  `warn_user`.

# Verification

`uv run --group test pytest tests/pp -q` plus `uvx prek run --all-files`.

# Gotchas

- `relative`'s all-zero-sample behavior is a deliberate departure from
  phyloseq, stated in the docstring `Notes` - do not "fix" it to return NaN.
- `tax_glom`'s output column order follows the *archetypes'* original
  position in `var_names`, not the order groups are first encountered while
  scanning lineages; `bt.datasets.toy()`'s taxa are contiguous per phylum, so
  `tests/pp/test_glom.py` builds its own interleaved fixture to cover this.
