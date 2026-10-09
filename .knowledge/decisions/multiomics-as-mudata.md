---
type: Decision
title: Multi-omics data is one MuData with fixed modality names
description: Several data types over the same samples are one MuData whose modalities are named taxa, function, function_by_taxon, metabolites and host; io.to_mudata keeps only the samples every modality has; a TreeData modality loses its tree in h5mu.
tags: [io, mudata, multiomics]
status: draft
paths: ["src/biotapy/io/_mudata.py"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T12:00:00Z }
commit: d23cd79
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Context
The spec stores several data types as MuData modalities `taxa`, `function`,
`metabolites` and `host`.[^spec] Phase 2 already made a function table a
MuData of two modalities, `function` and `function_by_taxon`
([function-tables-as-mudata](/decisions/function-tables-as-mudata.md)), and
MuData cannot nest one MuData in another as a modality. MuData also accepts
modalities whose samples only partly overlap: its global `obs` is the union and
`obsmap` marks the gaps, which no paired method (mmvec, a joint model) can use.

# Decision
- Modality names: `taxa` (a TreeData or AnnData of taxa), `function` and
  `function_by_taxon` side by side (a function table's two modalities,
  unchanged), `metabolites` and `host`. Functions that read a modality take
  its name as a keyword defaulting to these. Other names are allowed;
  nothing checks them.
- `io.to_mudata(modalities)` takes a mapping of name -> AnnData, keeps the
  samples every modality has in the first modality's order, copies each
  modality and warns once naming how many samples each loses. A MuData value
  raises `TypeError`: a function table goes in as `{**table.mod, ...}`.
- The global `obs` stays as MuData builds it (no columns); `mdata.pull_obs()`
  gathers them.
- Saving: `write_h5mu` drops a TreeData modality's `vart` and reads it back
  as AnnData ([tree-access](/contracts/tree-access.md), Gotchas). biotapy
  adds no writer; the docstring and the guide say to save that modality with
  `write_h5td` too.

# Rejected
- **Nested MuData** (`function` holding the function table): MuData's
  modalities must be AnnData.
- **Union of samples** (MuData's default): a sample missing from one modality
  becomes NaN rows that a paired method would have to drop anyway.
- **A biotapy `write_h5mu`** that writes trees beside the file: a second file
  format for one gotcha; reconsider when treedata or mudata supports it.
- **Rejecting unknown modality names**: blocks a `proteins` or `viruses`
  modality for no benefit.

# Consequences
- Every function that pairs modalities checks that they hold the same samples
  in the same order and names `io.to_mudata` as the fix.
- `datasets.hmp2` already follows the names (`function`,
  `function_by_taxon`, `taxa`).

[^spec]: Python Microbiome Toolkit development report, section Data model
