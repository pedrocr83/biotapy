---
type: Phase
title: Phase 4 - ML-ready and multi-omics (0.4)
description: MuData modality conventions, mmvec wrapper, leak-free scikit-learn transformers, an embedding plugin interface for microbiome foundation models, and a PyTorch loader.
tags: [roadmap, ml, multiomics]
status: stable
release: "0.4"
phase_state: in-progress
effort: ~4-6 weeks part-time
depends_on: [/roadmap/phase-3-stats.md]
paths: ["src/biotapy/ml/**", "src/biotapy/tl/**", "src/biotapy/io/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-08T17:09:16Z }
commit: bb393d4
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Goal
Differentiator 3: ML-ready by default, and multi-omics as MuData.[^spec]

# Entry criteria
- Phase 3 exit gate passed and 0.3 released. Nothing in this phase starts earlier.

# Design notes
- **Classes are allowed here only** because scikit-learn and PyTorch define
  class protocols (`BaseEstimator`/`TransformerMixin`, `torch.utils.data.Dataset`).
  One level of inheritance from the library base, no biotapy base classes.
- **Leakage comes from stateful steps**: prevalence/abundance filtering and
  any learned scaling. Per-sample transforms (relative, CLR) are stateless and
  only wrapped for pipeline convenience. The docs say which is which.
- **Embedding plugins** are discovered through the entry-point group
  `biotapy.embeddings`; each plugin exposes
  `embed(adata: AnnData, *, batch_size: int) -> np.ndarray` and results go to
  `obsm["X_<plugin>"]`. Model weights are never bundled.
- **Saving loses trees.** `write_h5mu` drops a TreeData modality's `vart` and
  reads back plain AnnData ([tree-access](/contracts/tree-access.md)). Task 4.1
  either writes tree-bearing modalities separately with `write_h5td`, or files
  an upstream issue with mudata/treedata and documents the limitation.

# Tasks
Expand into TDD steps (superpowers:writing-plans) when the phase starts.

- [ ] **4.1 MuData conventions** - modality names `taxa`, `function`,
  `metabolites`, `host`; `io.to_mudata(**modalities) -> MuData` aligning samples; decision concept.
  A Phase 2 function table is already a two-modality MuData (`function`,
  `function_by_taxon`), so `function` here cannot be a nested MuData: use
  `function` and `function_by_taxon` side by side
  ([function-tables-as-mudata](/decisions/function-tables-as-mudata.md)).
- [ ] **4.2 `tl.mmvec(mdata, *, microbes="taxa", metabolites="metabolites", seed=None) -> pd.DataFrame`** - wraps scikit-bio's mmvec.
- [ ] **4.3 `ml.PrevalenceFilter`, `ml.RelativeAbundance`, `ml.CLR`** - sklearn transformers over arrays (`fit` learns kept features); `sklearn.utils.estimator_checks.check_estimator` passes.
- [ ] **4.4 `ml.embed(adata, model, *, batch_size=64, inplace=False)`** - plugin loader over the entry-point group; one reference plugin (MGM or BiomeGPT, whichever has usable public weights).
- [ ] **4.5 `ml.to_torch(adata, *, label_key=None, layer=None) -> torch.utils.data.Dataset`** - sparse rows stay sparse until batch time; extra `torch`.
- [ ] **4.6 Leak-free CV notebook** - `Pipeline(PrevalenceFilter, CLR, classifier)` inside `cross_validate`, next to the leaky version with the difference shown.
- [ ] **4.7 Knowledge** - `Module` concept for `ml`; decision concept for the plugin interface.

# Exit gate
- [ ] Leak-free CV example executed in docs.
- [ ] One foundation model plugged in end to end.

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
