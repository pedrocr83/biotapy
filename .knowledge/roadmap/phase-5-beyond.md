---
type: Phase
title: Phase 5 - Beyond 0.4 (backlog)
description: User- and issue-driven backlog - time series, community state types, single-cell and spatial bridge, MCP server, federated meta-analysis. Not planned in detail.
tags: [roadmap, backlog]
status: draft
release: "0.5+"
phase_state: in-progress
depends_on: [/roadmap/phase-4-ml-multiomics.md]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T20:40:40Z }
commit: 5cdd6ee
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Rule
Nothing here gets tasks until 0.4 ships and a GitHub issue or user request
asks for it.[^spec] Each item then gets its own phase concept.

# Candidates
| Item | Notes |
|---|---|
| Time series (miaTime-style stability, trajectories) | needs a sample-time convention in `obs` first |
| Dirichlet multinomial mixtures, community state types | may belong in `tl` in 0.1-0.3 if users ask earlier |
| Single-cell / spatial bridge | scverse neighbours; revisit when a concrete use case exists |
| MCP server | exposes `tl`/`da` to agents; separate package likely |
| pyloseq importer | spec's mitigation if pyloseq gains traction; build when a pyloseq user asks |
| Federated meta-analysis | **blocked**: check employment IP and side-project clauses before any public work |

[^spec]: Python Microbiome Toolkit development report, sections Roadmap and Risks
