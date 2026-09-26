---
type: Decision
title: Knowledge in OKF, user docs in Sphinx
description: Contributor and agent knowledge (decisions, contracts, playbooks, roadmap) is an OKF v0.2 bundle in .knowledge/; user docs are a Sphinx site in docs/.
tags: [docs, okf]
status: stable
verified: { by: human:pedrocr83, at: 2026-09-26T09:40:17Z }
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: okf-spec
    resource: https://github.com/GoogleCloudPlatform/open-knowledge-format/blob/main/SPEC.md
    title: Open Knowledge Format v0.2 specification
    author: team:google-cloud
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Context
The spec wants minimal code comments with explanation moved to docs, ADRs, and
agent instructions that repeat the engineering rules.[^spec] It leaves open
Sphinx vs MkDocs. OKF gives agent-readable, diffable, progressively
disclosed knowledge with provenance and trust fields.[^okf-spec] The OKF copy in
`GoogleCloudPlatform/knowledge-catalog/okf` is a frozen snapshot; the canonical
home is `GoogleCloudPlatform/open-knowledge-format`.

# Decision
| Audience | Location | Format | Content |
|---|---|---|---|
| Contributors, coding agents | `.knowledge/` | OKF v0.2 | decisions (ADRs), contracts, playbooks, module maps, roadmap phases |
| Users | `docs/` | Sphinx + myst-nb (cookiecutter-scverse default) | getting started, user guide, method notes, tutorials, API reference, Coming from R |
| Everyone, per function | docstring | NumPy style | summary, params, returns, R equivalent, example, link to a `docs/` page |

- ADRs live only in `.knowledge/decisions/`; the docs site links to them rather
  than duplicating them.
- Maintenance rules: [maintain-knowledge](/playbooks/maintain-knowledge.md).
- OKF conformance is a pytest test, so CI enforces it.

# Rejected
- **MkDocs Material**: fine tool, but cookiecutter-scverse ships Sphinx and
  scverse ecosystem listing is easier on the standard stack.
- **ADRs in `docs/`**: agents would load the whole docs tree to find them;
  OKF `index.md` files give progressive disclosure.
- **Everything in OKF**: user docs need executed notebooks and autodoc, which
  Sphinx provides.

# Consequences
- Two doc roots; the table above is the only routing rule, repeated in `rules.md`.

[^okf-spec]: Open Knowledge Format v0.2 specification
[^spec]: Python Microbiome Toolkit development report, sections Documentation and Tooling
