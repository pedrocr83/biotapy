---
type: Contract
title: Tree access
description: Only biotapy/_core/_tree.py imports treedata or networkx, so a TreeData API change touches one file.
tags: [data-model, tree, dependencies]
status: stable
paths: ["src/biotapy/_core/_tree.py", "pyproject.toml"]
generated: { by: claude-code/claude-sonnet-5, at: 2026-10-03T00:50:26Z }
commit: 1ad037b
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
   edge list (`tree_from_edges`), reading the phylogeny (`get_tree`), Newick
   parsing (`tree_from_newick`), building a tree from an ape `phylo` edge
   matrix (`tree_from_phylo`), listing tips (`tree_tips`), relabeling
   tips (`relabel_tips`), and converting the phylogeny to a scikit-bio
   `TreeNode` (`get_skbio_tree`).
3. Newick parsing reuses scikit-bio (`TreeNode.read([text], convert_underscores=False)`);
   biotapy never writes its own parser.
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
- `tree_from_newick` always passes `convert_underscores=False`: scikit-bio's
  default turns unescaped `ASV_1` into `ASV 1`, corrupting ids.
- Malformed Newick surfaces as two scikit-bio exceptions, not one: with the
  default format sniffing, `TreeNode.read` raises `UnrecognizedFormatError`
  for most malformed text (the sniffer rejects it) and `NewickFormatError`
  only for text that sniffs as Newick but fails to parse (e.g. `(a:x,b:1);`).
  `tree_from_newick` catches exactly those two and raises `ValueError` naming
  its `argument` (readers pass `tree='<path>'`). Passing `format="newick"` to
  skip the sniffer is worse: some malformed text then leaks a bare
  `IndexError`, and a `FormatIdentificationWarning` is emitted.
- Internal-node labels (support values such as `0.95`) cannot survive as graph
  node names because they repeat across the tree; `tree_from_newick` drops them
  and assigns fresh `n0, n1, ...` names in preorder, skipping any value already
  used as a tip name.
- A missing branch length (no `:length` in the Newick text) becomes `nan`, not
  `0.0` or `None`, so downstream sum-of-branch-length code must handle NaN.
- `relabel_tips` raises `ValueError` when a new name already names a node
  that is not itself renamed: `nx.relabel_nodes` would silently merge the two
  nodes (a mixed `ASV1`/sequence DADA2 tree lost an edge this way). Swaps are
  allowed. `read_dada2` also requires every tip to be a DNA sequence.
- `make_treedata`'s tree/table alignment issues at most one `UserWarning` per
  call, even when both extra features and extra tips exist; it names both
  counts and keeps only the shared features, pruning the tree to kept tips
  plus ancestors. No shared feature is a hard `ValueError`, not a warning.
- `tree_from_phylo` reads ape's `phylo` edge-matrix convention directly:
  1-based node ids, tips numbered `1..len(tips)`, everything above that an
  internal node. It reuses `tree_from_newick`'s collision-free `n<i>` naming
  for internal nodes and the same `_require_unique_names` tip check, so a
  phyloseq-derived tree (`io.read_phyloseq`) and a Newick-derived one
  (`io.read_biom`/`read_qiime2`/`read_dada2`) name their internal nodes the
  same way. A missing `edge.length` becomes `nan`, same as `tree_from_newick`.
- scikit-bio's Faith PD and UniFrac need a root with at most two children, so
  `get_skbio_tree` splits a wider root with a zero-length node, and phyloseq
  instead roots an unrooted tree at a random tip.

[^treedata]: treedata 0.3.1 on PyPI
[^spec]: Python Microbiome Toolkit development report, section Risks
