---
type: Decision
title: What "methods agree" means in da.consensus
description: A feature is a consensus hit when at least min_methods methods call it at q < alpha and every calling method gives it the same sign; untested is not "not significant", opposite calls are a conflict, and the user picks the methods.
tags: [da, statistics, api]
status: stable
paths: ["src/biotapy/da/_consensus.py", "src/biotapy/da/_schema.py"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T19:56:32Z }
commit: 1153e6a
sources:
  - id: nearing
    resource: https://www.nature.com/articles/s41467-022-28034-z
    title: Nearing et al. 2022, Microbiome differential abundance methods produce different results across 38 datasets, Nature Communications
  - id: pelto
    resource: https://academic.oup.com/bib/article/26/2/bbaf130/8093585
    title: Pelto et al. 2025, Elementary methods provide more replicable results in microbial differential abundance analysis, Briefings in Bioinformatics
  - id: oma
    resource: https://microbiome.github.io/OMA/docs/devel/pages/differential_abundance.html
    title: Orchestrating Microbiome Analysis, Differential abundance chapter
---

# Context
Differential abundance methods disagree: on 38 datasets, 14 methods called very
different sets of features.[^nearing] Phase 3 runs several methods behind one
result schema (`contracts/data-model-slots`, DA results) and needs one rule for
"where they agree" that every reader of `da.consensus`'s table can check by hand.
The literature is split on the remedy: Nearing et al. recommend a consensus of
several methods; Pelto et al. and the OMA book recommend one elementary method
and warn that trying methods until one agrees is selective reporting.[^pelto][^oma]

# Decision
Per feature `f` and method `m`, over result tables the user computed:

- `called(f, m)` is `qvalue < alpha`, strictly: `q == alpha` is not called.
  Calls are recomputed from `qvalue`, never taken from a method's own flag
  (LinDA's `reject` uses `<=`; ANCOM-BC2's `Signif` sits on Holm by default).
- Every `qvalue` is Benjamini-Hochberg over the features that method tested,
  so "called" means the same false discovery rate in every column.
- `n_tested(f)` counts methods with a finite `pvalue`; a feature a method did
  not test (NaN, or absent from its table) is "not tested", never "not
  significant". `n_significant(f)` counts the calls.
- `consensus(f)` is true when `n_significant >= min_methods` and every calling
  method has the same non-zero `direction`. `conflict(f)` is true when calling
  methods disagree in sign; a conflict is never a consensus. `direction(f)` is
  the shared sign, 0 when there is no call or a conflict.
- `da.consensus(results, *, alpha=0.05, min_methods=2)` only combines tables:
  it never runs a method, so the list of methods is written in the user's code
  before any result is seen. Tables must come from different methods and
  compare the same `contrast`.

# Rejected
- **Intersection of all methods** (`min_methods = len(results)` as the only
  rule): one method that cannot test a feature (ANCOM-BC2 on a feature absent
  from one group) would veto it. Still available as `min_methods=len(results)`.
- **k of n without direction**: two methods calling a feature in opposite
  directions would count as agreement.
- **Agreement among all methods with an estimate** (every method's sign, called
  or not): a method that tested a feature but found nothing would block a
  consensus on a sign it never claimed.
- **Rank aggregation** (combining p-values or ranks): needs p-values that are
  comparable across methods, which they are not, and options no use case asks
  for (rules.md R2.3).
- **A one-call `da.consensus(adata, formula, *, methods, ...)`** that runs the
  methods: eight arguments (rules.md R5 allows six), a policy for a method that
  fails (R7.4 forbids dropping it silently), per-method options, and it hides
  which methods were chosen.

# Consequences
- The table is a robustness report, not a way to choose a method; the guide and
  `da.consensus`'s Notes say to fix the methods before looking.
- Methods that share a model (ANCOM-BC and ANCOM-BC2) agree more often for that
  reason alone; consensus counts methods, not independent evidence.
- Each method's call is stored as `significant_<method>`, so a reader of the
  table, or a plot of it, needs no `alpha`.

[^nearing]: Nearing et al. 2022, Microbiome differential abundance methods produce different results across 38 datasets, Nature Communications
[^pelto]: Pelto et al. 2025, Elementary methods provide more replicable results in microbial differential abundance analysis, Briefings in Bioinformatics
[^oma]: Orchestrating Microbiome Analysis, Differential abundance chapter
