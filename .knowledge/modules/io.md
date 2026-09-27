---
type: Module
title: io
description: File readers and writer for BIOM, QIIME 2 artifacts and DADA2 sequence tables, building every TreeData through _core.make_treedata.
resource: /src/biotapy/io/
paths: ["src/biotapy/io/**"]
tags: [io]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-27T07:28:41Z }
commit: 43d6efb
status: stable
---

# Responsibility

Owns the `bt.io.*` verbs that move data between files and a TreeData:
`read_biom`/`write_biom` (`_biom.py`), `read_qiime2` (`_qiime2.py`) and
`read_dada2` (`_dada2.py`), plus the private checked-join helper shared by the
last two (`_join.py`). Does NOT own phyloseq `.rds`/`.RData` reading
(`read_phyloseq`, Task 1.10, stage 2, not yet written) or downloaded example
datasets (`datasets.global_patterns`/`enterotype`, Task 1.11).

# Entry points

- `_biom.py:read_biom` / `_biom.py:write_biom` - BIOM 1.0 JSON and BIOM 2.1
  HDF5, both directions, through biom-format.
- `_biom.py:_biom_parts` - private; turns a loaded `biom.Table` into `(X, obs,
  var)`. Used by `read_biom` and reused by `read_qiime2` (see Gotchas).
- `_qiime2.py:read_qiime2` - a `.qza` feature table plus optional taxonomy,
  tree and QIIME 2 metadata TSV, with no QIIME 2 install.
- `_dada2.py:read_dada2` - a DADA2 `seqtab`/`seqtab.nochim` CSV or TSV, plus
  optional `assignTaxonomy`/`addSpecies` taxonomy and a Newick tree.
- `_join.py:_join_to` - private; the one checked reindex-join for a side file
  (taxonomy, metadata, taxa) onto a table's ids, used by `read_qiime2`
  (taxonomy, metadata) and `read_dada2` (taxa).

# Invariants

- Every reader builds its result only through `_core.make_treedata` and sets
  `x_kind` by calling `_core.infer_x_kind` on the freshly parsed matrix,
  never a literal. `_biom.py:read_biom`, `_qiime2.py:read_qiime2`,
  `_dada2.py:read_dada2`.
- biom-format stores observations (features) x samples; `_biom_parts`
  transposes exactly once so biotapy's samples-as-rows convention holds from
  there on. `_biom.py:_biom_parts`.
- `write_biom` writes every rank as a prefixed value (`k__`..`g__`, including
  a bare prefix for a missing rank) because BIOM's HDF5 reader drops empty
  taxonomy list entries on read; a bare prefix stays non-empty so position and
  rank survive the round trip. `_biom.py:_taxonomy_metadata`.
- BIOM has no null: `write_biom` writes a missing `obs` value as `""`
  (`_biom.py:_sample_metadata`), and `read_biom` maps `""` back to NaN on the
  way in (`_biom.py:_metadata_frame`), so a missing value survives a
  `write_biom` -> `read_biom` round trip as NaN, not `""`.
- `read_qiime2` recognizes a `.qza`'s payload by the path inside the zip
  (`<uuid>/data/<filename>`), not by parsing `metadata.yaml` - avoids a YAML
  runtime dependency for a check the payload path already gives for free.
  `_qiime2.py:_payload`.
- QIIME 2 metadata ID-header matching is mixed-case on purpose: modern headers
  (`sample-id`, `id`, ...) match case-insensitively, legacy headers
  (`#SampleID`, ...) must match exactly, mirroring QIIME 2's own metadata
  parser. `_qiime2.py:_is_id_header`, `_qiime2.py:_ID_HEADERS_ANY_CASE`,
  `_qiime2.py:_ID_HEADERS_EXACT`.
- `_join.py:_join_to` is the one checked join for a side file onto a table's
  ids: a duplicate id in the side file or zero shared ids raise `ValueError`
  naming the argument; some-but-not-all table ids covered gives one
  `UserWarning` naming the count, and ids in the side file but not the table
  are silently ignored (metadata commonly covers more samples than one
  table). `_join.py:_join_to`.
- `read_dada2` requires a supplied tree's tips to be DNA sequences (checked
  against `_dada2.py:_SEQUENCE`), never bare ASV ids, before relabeling them
  to `ASV<n>`: an ASV-id tip raises naming `tree=`, because biotapy's ASV
  numbering (table column order) has no way to be verified against a
  user-supplied id. `_dada2.py:_require_sequence_tips`, `_dada2.py:read_dada2`.
- `_dada2.py:_read_table` reads the seqtab's row-name column as text verbatim
  (`dtype={0: str}`, `keep_default_na=False`, `na_values=[""]`), so a
  numeric-looking or literal `"NA"` sample name survives unparsed; rank cells
  spelling `"NA"` are still normalized to NaN downstream by
  `_core.normalize_ranks`, not by `_read_table` itself.

# Dependencies

- [core](/modules/core.md): `make_treedata`, `infer_x_kind`, `split_lineage`,
  `normalize_ranks`, `tree_from_newick`, `tree_tips`, `relabel_tips`,
  `warn_user`, `TreeData`, `RANKS`, `as_csr`.
- Third-party: `biom-format` (`_biom.py`, and `_qiime2.py` via
  `_biom_parts`). No CPython 3.14 wheels yet (biocore/biom-format#1004): on
  3.14 it builds from source and needs a C compiler, noted in `README.md`.

# Verification

`uv run --group test pytest tests/io -q` plus `uvx prek run --all-files`.

# Gotchas

- `io/_qiime2.py` imports the private `_biom_parts` straight from its sibling
  `io/_biom.py` (`_qiime2.py:15`) rather than duplicating the BIOM-parsing
  logic. This is a deliberate exception to
  [module-boundaries](/contracts/module-boundaries.md) rule 1, now amended to
  allow a private import between topic files of the *same* subpackage; it
  would still be a violation across subpackages.
- A `#q2:types` row is padded through the same `_padded_row` logic as data
  rows and zipped with `strict=True` (`_qiime2.py:_metadata_rows`): a types
  row shorter than the header used to `zip(..., strict=False)` and silently
  drop the trailing declared-but-untyped columns from the result (fixed at
  Checkpoint B1, C2).
- A metadata column declared `numeric` (`#q2:types`) raises `ValueError`
  naming the column and the bad values on a non-numeric entry; only an
  *undeclared* column falls back to text when it fails to parse as numbers
  everywhere - `pd.to_numeric(errors="coerce")` alone would have silently
  turned bad values into NaN instead of rejecting the file the way QIIME 2
  itself does. `_qiime2.py:_typed`.
- biom-format ships no `py.typed`; rather than a blanket
  `disallow_untyped_calls = false`, mypy strict mode is narrowed with
  `untyped_calls_exclude = ["biom"]` (`[tool.mypy]` in `pyproject.toml`), so
  only calls into biom-format itself are exempt, not the rest of `io`.
