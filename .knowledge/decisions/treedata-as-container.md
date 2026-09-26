---
type: Decision
title: TreeData and MuData are the only containers
description: biotapy owns no data class; one data type is a TreeData, several are a MuData, every feature is a function over them.
tags: [architecture, data-model]
status: stable
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Context
A microbiome library needs counts, taxonomy, sample metadata, a phylogeny and
derived results in one object. phyloseq (R) and pyloseq (Python) each define
their own container. scverse already has AnnData, TreeData (AnnData plus
trees) and MuData (several AnnData modalities).[^spec]

# Decision
- Single data type: `treedata.TreeData`. Phylogeny lives in `vart["phylo"]`.
- Several data types (taxa, function, metabolites, host): `mudata.MuData`.
- biotapy defines no container class and no subclass of either.
- Slot layout is fixed by [data-model-slots](/contracts/data-model-slots.md).
- All tree access goes through [tree-access](/contracts/tree-access.md), because
  TreeData is young and its API may move.

# Rejected
- **Custom `Phyloseq`-style class** (pyloseq's route): no scverse interop, no
  h5ad/zarr, users must learn a new object. It is the gap biotapy exists to fill.
- **Plain AnnData with the tree in `uns`**: trees would not be subset with
  features, and `uns` is not a typed slot.

# Consequences
- Easy: scanpy/muon users already know the object; h5ad/zarr I/O is free.
- Hard: biotapy depends on TreeData's release cadence. Mitigated by pinning a
  version range in `pyproject.toml` and isolating tree calls in one file.

[^spec]: Python Microbiome Toolkit development report, sections Architecture and Data model
