---
type: Review
title: Spec review of plan.md
description: Gaps, contradictions and risks found in the development report, each with where it is resolved in this bundle.
tags: [roadmap, review]
status: stable
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Verdict
The spec is strong on positioning, engineering discipline and the performance
ladder. Its weak spots are the ones that bite during implementation: slot
semantics after transformations, what "matches R" means for stochastic
methods, and a few undefined terms. None blocks Phase 0.[^spec]

# Findings

| # | Finding | Severity | Resolution |
|---|---|---|---|
| 1 | Package name open, but repo is already `biotapy`; all five candidate names are free on PyPI (2026-09-26). | blocker for 0.0.1 | [package-name-biotapy](/decisions/package-name-biotapy.md), Phase 0 task 0.1 |
| 2 | "Pure by default" vs results stored in `obsm`/`obsp`: no rule for how `tl` returns or stores results. | high | [pure-by-default](/decisions/pure-by-default.md) |
| 3 | `tax_glom` example groups by rank value (phyloseq groups by full lineage), picks no representative, drops the tree and `uns`, and uses undefined `sum_by`/`first_by`; no rule for which slots survive a feature-changing op. | high | [data-model-slots](/contracts/data-model-slots.md) propagation and aggregation rules; Phase 1 task 1.5 |
| 4 | Golden tests cannot match R for rarefaction, NMDS, or permutation p-values (different RNGs); PCoA axes have arbitrary sign. | high | [r-golden-parity](/contracts/r-golden-parity.md) comparison table |
| 5 | `X` is "counts", but MetaPhlAn gives relative abundance and HUMAnN gives RPK/CPM; count-only methods (rarefy) would silently misbehave. | high | `uns["biotapy"]["x_kind"]` in [data-model-slots](/contracts/data-model-slots.md) |
| 6 | HUMAnN pathway abundances are not additive over taxa, so one table cannot hold both stratified and unstratified data; MuData is needed in 0.2, not 0.4. | medium | [phase-2-function](/roadmap/phase-2-function.md) design notes |
| 7 | KO -> pathway is many-to-many; a group-by sum is not enough. | medium | Phase 2 task 2.1 membership matrix |
| 8 | "Where methods agree" (DA consensus) is undefined. | medium | [phase-3-stats](/roadmap/phase-3-stats.md) design notes |
| 9 | "Private helpers only when shared by two or more functions" conflicts with the 30-line budget (splitting creates single-use helpers). | low | rules.md R4.3-R4.4 |
| 10 | "Reproduces the phyloseq tutorial" names no tutorial; the obvious one (`phyloseq-analysis.Rmd`) also uses CCA, DPCoA, `plot_tree` and `plot_net`, which are outside the 0.1 scope. "Top 30 phyloseq functions" is not listed. | medium | [phase-1-core](/roadmap/phase-1-core.md) exit gate: in-scope sections only, 31 functions listed |
| 11 | NMDS needs a non-metric MDS; scikit-bio has none, so scikit-learn becomes a core dependency in 0.1 (spec lists it only for `ml`). | medium | Phase 1 task 1.17, needs dependency approval |
| 12 | 0.1 scope (4 readers, ~17 functions, 4 plots, R golden infra, docs) in 6-8 part-time weeks is optimistic. | medium | Phase 1 split into slices 1A-1D with checkpoints |
| 13 | Leakage in `ml` comes from stateful steps (prevalence filtering), not from per-sample CLR; spec frames all transforms as leak risks. | low | [phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md) design notes |
| 14 | OKF lives at `GoogleCloudPlatform/open-knowledge-format`; the `knowledge-catalog/okf` copy is frozen. | info | [docs-okf-and-sphinx](/decisions/docs-okf-and-sphinx.md) cites the canonical repo |
| 15 | Spec says "CLAUDE.md or AGENTS.md"; only CLAUDE.md is created (imports rules.md). Other agents need an AGENTS.md pointer. | low | open: add `AGENTS.md` if another agent is used |
| 16 | CI matrix says Python 3.11-3.14, but treedata 0.3.1, mudata 0.4.1 and cookiecutter-scverse all require >= 3.12. | high | [phase-0-foundation](/roadmap/phase-0-foundation.md) global constraints |
| 17 | `unifrac` (Striped UniFrac) has no pip wheels (sdist needs conda headers, no Windows build); it cannot be an ordinary extra. | medium | [optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md): UniFrac via scikit-bio |
| 18 | scikit-bio 0.7's AnnData input reads `.X` with `np.asarray`, which fails on sparse matrices; wrappers must densify. | medium | rules.md R6.2; Phase 1 tasks 1.15-1.16 |
| 19 | MuData holds TreeData in memory but `write_h5mu` drops `vart`; multi-omics files silently lose the tree. | medium | [tree-access](/contracts/tree-access.md) gotchas; Phases 2 and 4 |
| 20 | Reading a phyloseq `.rds` via `rdata` is plausible (S4 supported) but unproven for a real phyloseq object. | medium | Phase 1 task 1.6 spike, then a decision |
| 21 | phyloseq defaults differ from the obvious Python calls: `rarefy_even_depth` samples with replacement, vegan `jaccard` is quantitative, weighted UniFrac is normalized. Naive golden tests would fail or, worse, pass for the wrong reason. | medium | Phase 1 tasks 1.14, 1.16 gotchas |

# Open questions still for the user
- Anchor cohort for tutorials (HMP2/IBDMDB proposed).
- scverse ecosystem listing: proposed after 0.1 ships.

[^spec]: Python Microbiome Toolkit development report
