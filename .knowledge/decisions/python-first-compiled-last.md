---
type: Decision
title: Python first, compiled code last
description: Pure Python/NumPy by default; delegate to compiled libraries; Numba then Rust only for a benchmarked hotspot; never new C/C++.
tags: [performance, architecture]
status: stable
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: skbio-numba
    resource: https://github.com/scikit-bio/scikit-bio/pull/2557
    title: scikit-bio Numba kernels PR
---

# Context
The heaviest microbiome math (UniFrac, PCoA, PERMANOVA) already runs as
compiled code in scikit-bio, unifrac and scikit-bio-binaries.[^spec] scikit-bio
reports 7x to 13x PERMANOVA speedups from its own Numba engine.[^skbio-numba]

# Decision
Escalation ladder. Each step only after the previous one fails a benchmark on
realistic data (5,000 samples x 50,000 features, sparse):

1. Vectorize with NumPy / SciPy sparse. Profile (py-spy, scalene) first.
2. Delegate to scikit-bio, unifrac, scikit-bio-binaries; expose their engine/GPU switch.
3. Numba kernel (optional dependency), selected by `engine="numba"`.
4. Rust via PyO3 + maturin, one crate in `rust/`, selected by `engine="rust"`.
5. New C/C++: never. C/C++ only through existing libraries.
6. GPU later via the array API (scikit-bio) and CuPy, not custom kernels.

Every compiled kernel obeys [engine-parity](/contracts/engine-parity.md).

# Rejected
- **Rust from day one**: raises the contributor bar and switches the build
  backend to maturin before any hotspot is proven.
- **New C/C++ extensions**: memory-unsafe and far harder to wheel than Rust.

# Consequences
- The build backend stays hatchling until a Rust crate exists.
- `asv` benchmarks are mandatory infrastructure from Phase 1, not an afterthought.

[^spec]: Python Microbiome Toolkit development report, section Performance
[^skbio-numba]: scikit-bio Numba kernels PR
