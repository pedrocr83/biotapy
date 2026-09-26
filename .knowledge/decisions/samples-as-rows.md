---
type: Decision
title: Samples are rows
description: Every matrix is samples x features (scverse orientation); importers transpose once, no other function checks orientation.
tags: [architecture, data-model, io]
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
phyloseq commonly stores `taxa_are_rows = TRUE`; BIOM and DADA2 outputs vary.
scverse (AnnData, scanpy, muon) always stores observations (samples) as rows.[^spec]

# Decision
- `X` is always samples x features.
- Only `biotapy.io` readers know about source orientation. Each reader
  transposes exactly once, before constructing the TreeData.
- No function outside `io` inspects or flips orientation.

# Rejected
- **Orientation flag on the object** (phyloseq's `taxa_are_rows`): every
  function would have to branch on it; it is the most common phyloseq bug source.

# Consequences
- Migration docs ("Coming from R") must state the flip prominently.
- R golden files are exported transposed so comparisons are direct; see
  [r-golden-parity](/contracts/r-golden-parity.md).

[^spec]: Python Microbiome Toolkit development report, section Data model
