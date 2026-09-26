---
type: Contract
title: Tree access
description: Only biotapy/_core/_tree.py imports treedata or networkx, so a TreeData API change touches one file.
tags: [data-model, tree, dependencies]
status: stable
paths: ["src/biotapy/_core/_tree.py", "pyproject.toml"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: b77a226
sources:
  - id: treedata
    resource: https://pypi.org/pypi/treedata/json
    title: treedata 0.3.1 on PyPI (2026-07-08, requires-python >=3.12)
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Statement
1. `src/biotapy/_core/_tree.py` is the only module that imports `treedata` or
   `networkx`. Everything else imports the `TreeData` type and tree helpers
   from `biotapy._core`.
2. It owns: constructing a TreeData (`make_treedata`), building a tree from an
   edge list (`tree_from_edges`), reading the phylogeny (`get_tree`), and, when
   first needed, Newick parsing and conversion to `skbio.TreeNode`.
3. Newick parsing reuses scikit-bio (`TreeNode.read([text])`); biotapy never
   writes its own parser.
4. `treedata` is pinned to `>=0.3.1,<0.4` in `pyproject.toml`.

# Why
TreeData is at 0.3.x and may change its API.[^treedata] The spec's mitigation
("wrap tree access in 2 or 3 helper functions") only works if nothing else
reaches past the helpers.[^spec]

# Enforced by
- ruff `TID251` banned-api for `treedata` and `networkx`, with a per-file
  ignore for `src/biotapy/_core/_tree.py` (Phase 0, task 0.4). Tests are exempt.

# Gotchas
- TreeData requires every tree to pass `nx.is_tree`, and with the default
  `alignment="leaves"` leaf names must be a subset of `var_names`.
- Subsetting prunes trees to kept leaves plus ancestors (unary nodes kept,
  branch lengths not merged).
- `MuData` holds TreeData modalities in memory, but `write_h5mu` drops `vart`
  and reading returns plain AnnData. Save tree-bearing modalities with
  `write_h5td` (matters from Phase 2).

[^treedata]: treedata 0.3.1 on PyPI
[^spec]: Python Microbiome Toolkit development report, section Risks
