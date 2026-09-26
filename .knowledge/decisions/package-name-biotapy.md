---
type: Decision
title: Package name biotapy
description: Distribution and import name is biotapy, hosted at github.com/pedrocr83/biotapy; confirm before the 0.0.1 placeholder release.
tags: [naming, release]
status: stable
verified: { by: human:pedrocr83, at: 2026-09-26T09:40:17Z }
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: pypi-check
    resource: https://pypi.org/pypi/biotapy/json
    title: PyPI JSON API lookup (HTTP 404 = name free)
    last_modified: 2026-09-26T08:20:00Z
---

# Context
The spec's shortlist is commensal, holobiont, taxonomia, microbiota. The repo
is already named `biotapy` with remote `git@github.com:pedrocr83/biotapy.git`.
On 2026-09-26 the PyPI JSON API returned 404 for all five names.[^pypi-check]

# Decision
- PyPI distribution: `biotapy`. Import: `import biotapy as bt`.
- GitHub: `pedrocr83/biotapy` (personal account) until an org is justified.
- Confirmed by the user on 2026-09-26 (Phase 0 task 0.1), with Python >= 3.12.

# Consequences
- Name is reserved only once 0.0.1 is published (Phase 0, task 0.8).

[^pypi-check]: PyPI JSON API lookup
