---
type: Decision
title: No bundled KEGG mapping files
description: Functional hierarchy mappings are downloaded on first use from user-provided or open sources (MetaCyc, eggNOG), never shipped in the wheel.
tags: [fn, licensing, datasets]
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
KEGG licensing restricts redistribution of its mapping files.[^spec]

# Decision
- `fn` hierarchy loaders take a path or URL; open sources are fetched and
  cached with pooch (same mechanism as `datasets`).
- No mapping file is committed to the repo or packaged in the wheel.
- Tests use a tiny synthetic hierarchy fixture, not real KEGG data.

# Consequences
- First use needs network access or a user-supplied file; docs say so.

[^spec]: Python Microbiome Toolkit development report, section Risks
