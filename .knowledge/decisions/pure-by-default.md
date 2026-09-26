---
type: Decision
title: Pure by default, one inplace convention for tl
description: io/pp return new objects and never mutate input; tl returns results, and inplace=True writes them to the documented slot; pl returns Axes.
tags: [api, conventions]
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
The spec says "pure by default; `inplace=True` only where scanpy users expect
it, and then consistently",[^spec] but also stores ordinations in `obsm` and
distances in `obsp`. Without a rule, every `tl` function would pick its own
behaviour. Confirmed by the user on 2026-09-26.

# Decision
| Module | Returns | Mutates input |
|---|---|---|
| `io`, `datasets` | new `TreeData` / `MuData` | n/a |
| `pp` | new `TreeData` | never |
| `tl` (default `inplace=False`) | the result (`pd.DataFrame`, `np.ndarray`, `pd.Series`) | never |
| `tl` with `inplace=True` | `None` | writes to the slot named in [data-model-slots](/contracts/data-model-slots.md), key overridable by `key_added` |
| `fn`, `da` | new object or result `pd.DataFrame` | never |
| `pl` | `matplotlib.axes.Axes` | never |

Every `tl` function supports both modes with identical semantics.

# Rejected
- **scanpy default (`copy=False`, mutate)**: contradicts the spec's purity rule.
- **Always return a copied TreeData from `tl`**: copies `X` on every call,
  which is hundreds of MB at 5,000 x 50,000.

# Consequences
- Tests assert the input object is unchanged after every `pp`/`tl` call.

[^spec]: Python Microbiome Toolkit development report, section Function-level implementation
