# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/
[semantic versioning]: https://semver.org/

## [Unreleased]

## [0.1.0] - 2026-10-03

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
- `bt.pp.filter_features` / `bt.pp.filter_samples`: keep features by
  prevalence and total reads, and samples by depth.
- `bt.pp.rarefy`: subsample every sample to the same depth, without
  replacement.
- `bt.datasets.esophagus`: phyloseq's esophagus dataset, with a tree.
- `bt.tl.alpha`: observed features, Shannon, Simpson, Chao1 and Faith's PD.
- `bt.tl.beta` / `bt.tl.unifrac`: Bray-Curtis and Jaccard distances, and
  unweighted and weighted UniFrac.
- `bt.tl.pcoa` / `bt.tl.nmds`: principal coordinates and non-metric
  multidimensional scaling of a stored distance matrix.
- `bt.tl.permanova`: PERMANOVA of a stored distance matrix.
- R golden parity tests for `pp.relative` and `pp.tax_glom` against phyloseq
  on real data.
- New runtime dependencies: `rdata`, `xarray`, `pooch`, `scikit-bio`,
  `biom-format`, `treedata`, `networkx`, `scipy`, `pandas`.
- New runtime dependency: `scikit-learn>=1.8`, for `tl.nmds`.
- `bt.pl.bar`, `bt.pl.heatmap`, `bt.pl.richness`, `bt.pl.ordination` and
  `bt.pl.scree`: phyloseq's plots, drawn from the results `tl` and `pp`
  store, returning matplotlib axes.
- A Coming-from-R page, generated from every public docstring's R
  equivalent, and a test that keeps those docstrings parseable.
- Tutorials: getting started, and phyloseq's analysis vignette redone with
  biotapy. Every notebook now runs on each docs build.
- asv benchmarks on a synthetic 5,000 x 50,000 table, with the 0.1
  baselines in the docs.
- New runtime dependency: `matplotlib>=3.8`, for `pl`, imported only when a
  plot is drawn.
- New runtime dependency: `threadpoolctl>=3.5`, already required by
  scikit-learn, to limit `tl.permanova`'s OpenMP threads.

### Changed

- Feature-changing functions (`pp.filter_features`, `pp.rarefy`,
  `pp.tax_glom`) keep `X` and drop the ordination summaries
  `uns["biotapy"]["pcoa"]` and `["nmds"]` along with the other derived slots.
- CI builds the docs, executing every notebook, and fails when total line
  coverage drops below 90%.
- `bt.tl.permanova` runs scikit-bio's OpenMP permutation test on one thread: on
  a machine with busy cores, its default thread pool was up to thousands of
  times slower.
- Package metadata declares the license as the SPDX expression `BSD-3-Clause`
  (PEP 639) and adds development status, audience and topic classifiers.
- The sdist no longer ships the tests that check the repository itself (CI
  config and knowledge bundle), so its test suite passes on its own.

## [0.0.1] - 2026-09-26

### Added

- Package skeleton, CI, documentation site and development rules. No public functions yet; this release reserves the name.
