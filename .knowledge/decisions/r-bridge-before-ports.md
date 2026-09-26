---
type: Decision
title: R bridge before native ports
description: R-only DA methods (MaAsLin 3, ANCOM-BC2, ALDEx2) ship first through an optional rpy2 bridge; native ports only for the most used, validated against R.
tags: [da, r, dependencies]
status: stable
verified: { by: human:pedrocr83, at: 2026-09-26T09:40:17Z }
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Context
Several differential abundance methods exist only in R. Porting them is large,
error-prone work, and the spec lists "port MaAsLin 3 natively or rpy2 only?"
as open.[^spec]

# Decision
- Phase 3 ships `da` bridges over rpy2 (extra `r`), each with the same output
  schema as native methods.
- scikit-bio's ANCOM-BC is the native default where it exists.
- A native port is considered only after 0.3, for a method with demonstrated
  demand, and must match the R bridge on the benchmark dataset.

# Rejected
- **Native MaAsLin 3 port in 0.3**: blows the 4-week Phase 3 budget.

# Consequences
- CI needs a job with R available for bridge tests; other jobs skip them via a
  pytest marker `r`.

[^spec]: Python Microbiome Toolkit development report, sections Risks and Open questions
