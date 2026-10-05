---
type: Phase
title: Phase 3 - Differential abundance consensus (0.3)
description: CLR and PhILR transforms; ANCOM-BC, LinDA, ALDEx2 and MaAsLin 3 behind one result schema; a consensus runner that reports where methods agree.
tags: [roadmap, da, pp]
status: stable
release: "0.3"
phase_state: in-progress
effort: ~4 weeks part-time
depends_on: [/roadmap/phase-2-function.md]
paths: ["src/biotapy/da/**", "src/biotapy/pp/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-05T13:25:58Z }
commit: 7df6215
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: maaslin3
    resource: https://www.nature.com/articles/s41592-025-02923-9
    title: MaAsLin 3, Nature Methods
---

# Goal
Differentiator 2: one call runs several DA methods and reports agreement.[^spec]

# Entry criteria
- Phase 2 exit gate passed and 0.2 released.
- [r-bridge-before-ports](/decisions/r-bridge-before-ports.md) confirmed by the user.

# Design notes (resolve before task 3.3)
- **"Agree" is undefined in the spec.** Proposed: a feature is a consensus hit
  when at least `min_methods` methods call it significant at `q < alpha` *and*
  all of those agree on effect direction. Written as a `Decision` concept before coding.
- **One result schema** for every method, a `pd.DataFrame` indexed by feature:
  `effect, se, pvalue, qvalue, direction, method`. Methods that lack `se`
  fill `NaN`, never drop the column.

# Tasks
Expand into TDD steps (superpowers:writing-plans) when the phase starts.

- [ ] **3.1 `pp.clr(adata, *, pseudocount=0.5, layer=None) -> AnnData`** - writes `layers["clr"]`;
  property: rows sum to 0. Thin wrapper over `skbio.stats.composition.clr`.
- [ ] **3.2 `pp.philr(tdata, *, pseudocount=0.5) -> AnnData`** - balances over the
  binary tree via `_core` tree helpers; Python over flattened node arrays first
  ([python-first-compiled-last](/decisions/python-first-compiled-last.md)); golden vs R `philr`.
- [ ] **3.3 Result schema** - `da/_schema.py`: `validate_result(df) -> pd.DataFrame`, used by every method; decision concept for "agree".
- [ ] **3.4 `da.ancombc(adata, formula, *, alpha=0.05) -> pd.DataFrame`** - wraps scikit-bio ANCOM-BC; golden vs R `ANCOMBC`.
- [ ] **3.5 `da.linda(adata, formula, *, alpha=0.05) -> pd.DataFrame`** - native (CLR + linear model + bias correction via statsmodels); golden vs R `LinDA`. statsmodels is a new dependency: ask.
- [ ] **3.6 `da.aldex2(adata, group, *, mc_samples=128, seed=None) -> pd.DataFrame`** - rpy2 bridge (extra `r`); marker `r`.
- [ ] **3.7 `da.maaslin3(adata, formula, *, seed=None) -> pd.DataFrame`** - rpy2 bridge (extra `r`); marker `r`.[^maaslin3]
- [ ] **3.8 `da.consensus(adata, formula, *, methods, alpha=0.05, min_methods=2, n_jobs=1) -> pd.DataFrame`** -
  joblib over methods; columns per method plus `n_significant`, `consensus`.
- [ ] **3.9 `pl.consensus(...)`** - UpSet-style or dot plot of agreement.
- [ ] **3.10 Docs** - method notes per DA method in `docs/methods/`; user guide page on consensus; `Module` concept for `da`.
- [ ] **3.11 CI** - job with R + Bioconductor for `-m r` tests.

# Exit gate
- [ ] Consensus report on one benchmark dataset, executed notebook in docs.
- [ ] Per-method agreement with the R reference, per [r-golden-parity](/contracts/r-golden-parity.md).

# Risks
- rpy2 on Windows/macOS CI is fragile -> R job runs on Linux only; bridges documented as Linux/macOS-first.

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
[^maaslin3]: MaAsLin 3, Nature Methods
