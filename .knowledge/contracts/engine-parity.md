---
type: Contract
title: Engine parity
description: A compiled kernel is a drop-in behind an existing public function via `engine=`; the Python engine is the test oracle and both must match within a stated tolerance.
tags: [performance, testing]
status: stable
paths: ["src/biotapy/**/*.py", "rust/**", "benchmarks/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T21:29:30Z }
commit: ebe0cbc
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Statement
1. An engine is selected by `engine: Literal["python", "numba", "rust"] = "python"`
   (or a delegated library's own switch, exposed under the same name).
2. Adding an engine never changes the public signature or return type.
3. The `"python"` engine is never removed; it is the reference oracle.
4. The same test function runs across every available engine
   (`@pytest.mark.parametrize("engine", available_engines())`) with an explicit
   `rtol`/`atol` written in the test.
5. Merge requires an `asv` benchmark on the 5,000 x 50,000 sparse synthetic
   dataset (`benchmarks/benchmarks/_data.py:synthetic`, baselines in
   `docs/performance.md`) showing >= 5x speedup, plus a note on the function's
   docs page.
6. An engine whose optional dependency is missing raises `ImportError` naming
   the extra; it never silently falls back to another engine.

# Why
Without an oracle, a fast kernel that is subtly wrong ships unnoticed; without
the 5x bar, compiled code accumulates for marginal gains and raises the
contributor bar for nothing. See [python-first-compiled-last](/decisions/python-first-compiled-last.md).

# Enforced by
- Parametrized engine tests (first instance: Phase 1 perf track, if any).
- PR template checkbox "benchmark attached" (Phase 0, task 0.6).
- Not otherwise automated: reviewers must check the asv result.
