---
type: Decision
title: No bundled KEGG mapping files
description: Function hierarchies come from the user's local files or from ENZYME (CC BY 4.0, `bt.datasets.enzyme`); KEGG and MetaCyc are never shipped or fetched.
tags: [fn, licensing, datasets]
status: stable
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-03T17:00:00Z }
commit: 8830258
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Context
KEGG licensing restricts redistribution of its mapping files.[^spec]

# Decision
- `fn.load_hierarchy` reads local files only, no URL. The one built-in
  download is ENZYME's EC hierarchy, in `datasets` with the other pooch
  fetches. No KEGG or MetaCyc mapping is committed, packaged, downloaded or
  used as a fixture: MetaCyc has been subscription-only since 2024, and KEGG's
  REST API is for academic users only. Test fixtures are synthetic maps and a
  CC BY 4.0 ENZYME excerpt.

# Consequences
- KEGG and MetaCyc users point `load_hierarchy` at files they are licensed to
  use.

[^spec]: Python Microbiome Toolkit development report, section Risks
