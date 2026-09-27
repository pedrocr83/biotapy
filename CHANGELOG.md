# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/
[semantic versioning]: https://semver.org/

## [Unreleased]

### Added

- `bt.io.read_biom` / `bt.io.write_biom`: read and write BIOM tables, both the
  1.0 (JSON) and 2.1 (HDF5) dialects.
- `bt.io.read_qiime2`: read QIIME 2 `.qza` artifacts (table, taxonomy, tree,
  sample metadata), with no QIIME 2 install needed.
- `bt.io.read_dada2`: read a DADA2 sequence table (CSV, TSV or `.rds`), with
  optional taxonomy and tree.
- `bt.io.read_phyloseq`: read a phyloseq object saved from R (`.rds`/`.RData`),
  with no R or rpy2 install needed.
- `bt.datasets.toy`: a tiny in-memory dataset for examples and tests.
- `bt.datasets.global_patterns` / `bt.datasets.enterotype`: phyloseq's example
  datasets, downloaded once and cached on disk with pooch.
- `bt.pp.relative`: per-sample relative abundance, added as a new layer.
- `bt.pp.tax_glom`: aggregate features that share a lineage to a taxonomic
  rank.
- R golden parity tests for `pp.relative` and `pp.tax_glom` against phyloseq
  on real data.
- New runtime dependencies: `rdata`, `xarray`, `pooch`, `scikit-bio`,
  `biom-format`, `treedata`, `networkx`, `scipy`, `pandas`.

## [0.0.1] - 2026-09-26

### Added

- Package skeleton, CI, documentation site and development rules. No public functions yet; this release reserves the name.
