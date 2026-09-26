---
type: Module
title: pp (preprocessing)
description: Pure transforms over AnnData/TreeData that scale abundances per sample or aggregate features along the taxonomy.
resource: /src/biotapy/pp/
paths: ["src/biotapy/pp/**"]
tags: [pp]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T14:08:28Z }
commit: 0fdbd4d
status: stable
---

# Responsibility

Owns the `bt.pp.*` verbs that transform an AnnData/TreeData's abundance
table: `relative` (adds a layer, keeps every feature) and `tax_glom`
(aggregates features to a taxonomic rank, drops derived slots). Does NOT own
filtering or rarefaction (`pp.filter_features`, `pp.filter_samples`,
`pp.rarefy` - Slice 1C, not yet written) or any diversity/ordination
computation (`tl`, later phases).

# Entry points

- `_transform.py:relative` - per-sample scaling to relative abundance;
  all-zero samples stay all-zero, which differs from phyloseq's
  `transform_sample_counts` (returns `NaN` there).
- `_glom.py:tax_glom` - aggregate features sharing a lineage down to `rank`
  into their most abundant member.

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

# Dependencies

- [core](/modules/core.md): `as_csr`, `sum_by`, `argmax_by`, `split_ranks`,
  `feature_subset`, `add_provenance`.

# Verification

`uv run --group test pytest tests/pp -q` plus `uvx prek run --all-files`.

# Gotchas

- `relative`'s all-zero-sample behavior is a deliberate departure from
  phyloseq, stated in the docstring `Notes` - do not "fix" it to return NaN.
- `tax_glom`'s output column order follows the *archetypes'* original
  position in `var_names`, not the order groups are first encountered while
  scanning lineages; `bt.datasets.toy()`'s taxa are contiguous per phylum, so
  `tests/pp/test_glom.py` builds its own interleaved fixture to cover this.
