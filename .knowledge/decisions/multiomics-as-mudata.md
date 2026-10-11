---
type: Decision
title: Multi-omics data is one MuData with fixed modality names
description: Several data types over the same samples are one MuData whose modalities are named taxa, function, function_by_taxon, metabolites and host; io.to_mudata keeps only the samples every modality has; `io.write_h5mu` and `io.read_h5mu` keep a TreeData modality's trees in h5mu, which plain mudata drops.
tags: [io, mudata, multiomics]
status: stable
paths: ["src/biotapy/io/_mudata.py", "src/biotapy/tl/_mmvec.py"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T23:26:30Z }
commit: 5123331
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
  its name as a keyword defaulting to these (`tl.mmvec(microbes="taxa",
  metabolites="metabolites")`). Other names are allowed; nothing checks them.
- `io.to_mudata(modalities)` takes a mapping of name -> AnnData, keeps the
  samples every modality has in the first modality's order, copies each
  modality and warns once naming how many samples each loses. A MuData value
  raises `TypeError`: a function table goes in as `{**table.mod, ...}`.
- The global `obs` stays as MuData builds it (no columns); `mdata.pull_obs()`
  gathers them.
- Saving: mudata's `write_h5mu` writes a TreeData modality as an AnnData, so
  its trees are lost ([tree-access](/contracts/tree-access.md), Gotchas).
  `io.write_h5mu(mdata, path)` / `io.read_h5mu(path)` keep them (task 5.1,
  user request 2026-10-11, reversing the earlier "no biotapy writer"): the
  file is written by mudata, then each TreeData modality's `obst`, `vart`,
  `label`, `allow_overlap` and `alignment` are added under its group
  (`_core/_tree.py:write_tree_slots`). Plain `mudata.read_h5mu` still opens the
  file and ignores them. h5mu only, not zarr. Upstream:
  [mudata#210](https://github.com/scverse/mudata/issues/210) (PR #211, a
  duck-typed `_write_mudata_extras` / `_read_mudata_extras` hook) and
  [treedata#102](https://github.com/YosefLab/treedata/issues/102) (PR #103).

# Rejected
- **Nested MuData** (`function` holding the function table): MuData's
  modalities must be AnnData.
- **Union of samples** (MuData's default): a sample missing from one modality
  becomes NaN rows that a paired method would have to drop anyway.
- **Trees beside the file** (a `.h5td` per modality): a second file to keep in
  step; the trees live in the one `.h5mu` instead.
- **Rejecting unknown modality names**: blocks a `proteins` or `viruses`
  modality for no benefit.

# Consequences
- Every function that pairs modalities checks that they hold the same samples
  in the same order and names `io.to_mudata` as the fix (`tl.mmvec`,
  `tl/_mmvec.py`).
- `datasets.hmp2` already follows the names (`function`,
  `function_by_taxon`, `taxa`).

[^spec]: Python Microbiome Toolkit development report, section Data model
