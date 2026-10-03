---
type: Decision
title: Function tables are a two-modality MuData
description: A HUMAnN-style function table is a MuData with a community modality and a stratified modality, adopted in Phase 2 instead of Phase 4; mudata is a runtime dependency.
tags: [fn, io, mudata, dependencies]
status: stable
paths: ["src/biotapy/_core/_function.py", "src/biotapy/io/_humann.py", "src/biotapy/io/_picrust2.py", "src/biotapy/fn/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-03T19:02:00Z }
commit: fbfeb99
sources:
  - id: research
    resource: ../roadmap/phase-2-function.md
    title: Phase 2 plan, design notes (function tables)
---

# Context
A pathway's community abundance is not the sum of its per-taxon strata, so the
stratified rows cannot be derived from the community rows and summing one table
that holds both double counts. Phase 2 needs both tables, aligned by sample,
from one `read_humann` call (`_core/_function.py:make_function_mudata`).
PICRUSt2's two files fill the same two modalities: `io.read_picrust2` takes
the unstratified table for `"function"` and the long contribution table
(`contrib=`) for `"function_by_taxon"` (`io/_picrust2.py:read_picrust2`), so
the stratified modality is empty without `contrib`. Its per-ASV trait table
describes genomes, not samples, so it is a DataFrame and not a modality
(`io/_picrust2.py:read_picrust2_traits`).

# Decision
`io.read_humann` and `io.read_picrust2` return a `MuData` with modalities `"function"` and
`"function_by_taxon"`, each an AnnData with its own copy of `obs`
(`_core/_function.py:FUNCTION_KEY`, `BY_TAXON_KEY`;
[data-model-slots](/contracts/data-model-slots.md), Function tables). Both
always exist, either may have zero features. `fn` verbs take a modality
(`func_glom`) or the MuData (`renorm`). `mudata>=0.4` is a runtime dependency
(`pyproject.toml`, `dependencies`), approved 2026-10-03
([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)). It
was brought forward from Phase 4, where multi-omics needs it anyway.

# Rejected
- **`dict[str, AnnData]`**: no dependency, but nothing keeps the two tables'
  samples aligned, there is no single-file IO (two h5ad files), and Phase 4
  wants MuData regardless.
- **Stratified table in `uns["biotapy"]`**: `uns` is not subset by AnnData
  slicing, so it goes silently stale after `pp.filter_samples`; it also breaks
  the contract that `uns["biotapy"]` holds metadata only.
- **`obsm` or `layers` slot**: layers need the same `var` as `X`, and the
  stratified rows have different features; `obsm` would need the (function,
  taxon) labels stored elsewhere, back in `uns`, with the same staleness.
- **One AnnData holding both kinds of row** (HUMAnN's own file layout): it
  subsets correctly but makes `X.sum()` double count, every function must
  filter on a flag, and `renorm`/`func_glom` would have to select rows first.
  The two-modality form makes the footgun unrepresentable.

# Consequences
- **h5mu drops a TreeData modality's tree**: `write_h5mu` then `read_h5mu`
  returns plain AnnData without `vart`
  ([tree-access](/contracts/tree-access.md), Gotchas). No function modality
  has a tree, so Phase 2 is unaffected.
- **mudata joins the dependencies**, with `scverse-misc[settings]` (pydantic-settings,
  python-dotenv) as new transitive packages; mypy treats it like other untyped
  libraries (`pyproject.toml`, `[[tool.mypy.overrides]]`).
- **Callers pass `mdata["function"]`, not the MuData**, to AnnData-typed verbs
  such as `fn.func_glom`; MuData does not guarantee modality
  alignment after a user edits one modality, so `fn.renorm` checks it
  (`fn/_renorm.py:_function_modalities`).
- **Forward note for Phase 4 task 4.1** ([phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md)):
  a function table is itself a two-modality MuData, so a multi-omics MuData
  needs `function` and `function_by_taxon` side by side with `taxa`,
  `metabolites` and `host`, not a nested `function` MuData. Task 4.1's modality
  names must not reuse `function` for anything else, and `io.to_mudata` must
  accept a function table's two modalities as two entries.
