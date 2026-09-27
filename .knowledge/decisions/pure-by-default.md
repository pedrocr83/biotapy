---
type: Decision
title: Pure by default, one inplace convention for tl
description: io/pp return new objects and never mutate input; tl returns results, and inplace=True writes them to the documented slot; pl returns Axes.
tags: [api, conventions]
status: stable
verified: { by: human:pedrocr83, at: 2026-09-27T19:50:20Z }
generated: { by: claude-code/claude-sonnet-5, at: 2026-09-27T18:29:04Z }
commit: 795dc19
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
| `tl` with `inplace=True` | `None` | writes to the slot named in [data-model-slots](/contracts/data-model-slots.md) |
| `fn`, `da` | new object or result `pd.DataFrame` | never |
| `pl` | `matplotlib.axes.Axes` | never |

Every `tl` function that returns per-sample or per-pair values supports both modes with
identical semantics; `tl.permanova` is the exception (see Consequences).

# Rejected
- **scanpy default (`copy=False`, mutate)**: contradicts the spec's purity rule.
- **Always return a copied TreeData from `tl`**: copies `X` on every call,
  which is hundreds of MB at 5,000 x 50,000.

# Consequences
- Tests assert the input object is unchanged after every `pp`/`tl` call.
- `tl.permanova` returns a test result, not per-sample or per-pair values, so
  it has no slot and no `inplace`. No `tl` function takes `key_added` until a
  use case needs one (rules.md R2.3).

[^spec]: Python Microbiome Toolkit development report, section Function-level implementation
