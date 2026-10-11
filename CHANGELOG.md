# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/
[semantic versioning]: https://semver.org/
[scikit-bio-2631]: https://github.com/scikit-bio/scikit-bio/issues/2631

## [Unreleased]

### Added

- `bt.io.write_h5mu` and `bt.io.read_h5mu`: save a `MuData` to an `.h5mu` file
  and read it back with each TreeData modality's trees kept (mudata alone
  writes a TreeData as an AnnData). biotapy now declares `h5py` as a
  dependency.

## [0.4.1] - 2026-10-10

### Fixed

- `bt.ml.embed(..., "mgm")`: the first embedding in a process no longer
  differs from later ones by up to 2e-4 on CPUs running several threads;
  biotapy fills oneMKL's CPU-type cache with one serial call before MGM's
  forward pass (pytorch/pytorch#188792).

## [0.4.0] - 2026-10-10

### Added

- `bt.io.to_mudata`: combine data types measured on the same samples into one
  `MuData`, keeping the samples every modality has (in the first one's order,
  with a warning naming how many each lost). Modality names: `taxa`,
  `function`, `function_by_taxon`, `metabolites`, `host`.
- `bt.tl.mmvec`: which metabolites go with which microbe (mmvec, Morton et
  al. 2019), fitted by scikit-bio, as a microbes x metabolites table of
  row-centred log probabilities.
- `bt.ml.PrevalenceFilter` and `bt.ml.CLR`: scikit-learn transformers, so a
  prevalence filter and CLR are fitted inside each cross-validation fold.
  They pass scikit-learn's estimator checks and give `bt.pp.filter_features`'s
  and `bt.pp.clr`'s results. Relative abundance is scikit-learn's own
  `Normalizer(norm="l1")`.
- `bt.ml.to_torch`: a PyTorch dataset over a table's samples, densifying one
  row at a time, with integer or float labels from an `obs` column. It needs
  the new extra `biotapy[torch]` (`torch>=2.9`).
- `bt.ml.embed`: one embedding per sample from a pretrained model, returned or
  stored in `obsm["X_<model>"]`. Models are plugins: a package registers one
  in the entry-point group `biotapy.embeddings`, and biotapy checks what it
  returns.
- MGM, the Microbial General Model (Zhang et al. 2026), as the first such
  model, behind the new extra `biotapy[mgm]` (`torch>=2.9`,
  `transformers>=5`). Its pretrained weights (MIT) are downloaded once from
  the `microformer-mgm` 0.5.8 wheel on PyPI and checked against its SHA-256;
  biotapy's embeddings are within 2e-6 of MGM's own code.
- `bt.datasets.biocrust`: mmvec's soil example, microbes and metabolites of
  a desert biocrust after wetting, downloaded and cached on first use.
- Guide pages for multi-omics and machine learning, and three tutorials:
  multi-omics with mmvec, leak-free cross-validation on the HMP2 cohort, and
  MGM embeddings of GlobalPatterns.

### Changed

- The pseudocount warning of `bt.pp.clr` and `bt.pp.philr` ends with the
  step that gave it, `(pp.clr)` or `(pp.philr)`.
- `bt.pp.filter_features` raises `TypeError` naming `min_prevalence` or
  `min_total` when the threshold is a bool or not a number.
- `bt.pp.tax_glom`, `bt.fn.func_glom` and `bt.pl.bar`'s `fill` sum a table
  of bools or integers narrower than 64 bits in 64 bits, as NumPy's `sum`
  does (unsigned stays unsigned): an `int8` table no longer wraps around, and
  a bool table is counted instead of saturating at `True`. An `int32` table
  now gives `int64` sums.
- `bt.da.ancombc2` documents that scikit-bio 0.7.4 can fail to estimate the
  bias with two samples in the reference level, where R returns results;
  biotapy raises naming scikit-bio
  ([scikit-bio#2631][scikit-bio-2631]).
- The Coming-from-R page marks the phyloseq calls biotapy does not cover yet
  "not in 0.4".

## [0.3.0] - 2026-10-07

### Added

- `bt.pp.clr`: centred log-ratio transform into `layers["clr"]`, equal to
  `vegan::decostand(x, "clr", pseudocount = 0.5)`.
- `bt.pp.philr`: PhILR balances along the phylogeny into `obsm["X_philr"]`,
  equal to `philr::philr` with its default weights.
- `bt.da.linda`, `bt.da.ancombc2`: differential abundance by LinDA (a port
  of `MicrobiomeStat::linda`) and ANCOM-BC2 (through scikit-bio), both in
  Python and compared with R in golden tests.
- `bt.da.aldex2`, `bt.da.maaslin3`: ALDEx2 and MaAsLin 3's abundance model,
  run in R through rpy2. They need R, the R package, and the new extra
  `biotapy[r]` (`rpy2>=3.6.8`, GPL-2.0-or-later, built against your R).
- One result table for every method: log2 `effect`, `se`, `pvalue`,
  Benjamini-Hochberg `qvalue` over the features the method tested,
  `direction`, `method` and `contrast`. Methods take `group`, `covariates`
  and `reference` instead of a formula, and never filter features.
- `bt.da.consensus` and `bt.pl.consensus`: count, per feature, the methods
  that call it and whether they agree on its direction, as a table and as a
  dot matrix.
- A differential abundance guide, a page per method, and a tutorial on the
  GlobalPatterns genera.
- asv benchmarks for `pp.philr`, `da.linda` and `da.ancombc2`, with their
  baselines in the docs.

### Changed

- The readers infer `x_kind` `"counts"` only for non-negative whole numbers,
  and the functions that need counts (`pp.rarefy`, `tl.alpha`'s
  `observed_features` and `chao1`, weighted `tl.unifrac`) refuse a table with
  a negative value.
- `bt.pl.bar` and `bt.pl.heatmap` raise `ValueError` for a table with a
  negative value, such as `layers["clr"]`, instead of drawing it as
  abundances.
- The Coming-from-R page marks the phyloseq calls biotapy does not cover yet
  "not in 0.3".

## [0.2.0] - 2026-10-05

### Added

- `bt.io.read_humann`: read HUMAnN 3 and 4 tables (gene families, reactions
  or pathway abundance, also after regrouping or renormalising) into a
  `MuData` with a community modality, `"function"`, and a per-taxon one,
  `"function_by_taxon"`.
- `bt.io.read_metaphlan`: read MetaPhlAn 3 and 4 profiles, one or merged, as
  relative abundances of the leaf clades with the seven taxonomic ranks.
- `bt.io.read_picrust2` / `bt.io.read_picrust2_traits`: read PICRUSt2's
  predicted metagenomes and pathways, with the per-ASV contributions as the
  per-taxon modality, and its per-ASV gene copy numbers. PICRUSt2 itself is
  not needed.
- `bt.fn.load_hierarchy`: read a function hierarchy (HUMAnN or PICRUSt2
  mapping files, or your own) from a local file.
- `bt.fn.func_glom`: aggregate functions to one level of a hierarchy, as
  `humann_regroup_table` does.
- `bt.fn.renorm`: renormalise a function table to copies per million or
  relative abundance, as `humann_renorm_table` does. `func_glom` and `renorm`
  are tested against HUMAnN 3.9's own output.
- `bt.fn.contributions` and `bt.pl.contributions`: one function's abundance
  per taxon in every sample, as a table and as stacked bars.
- `bt.fn.functional_redundancy`: taxonomic diversity, functional diversity
  and functional redundancy per sample (Tian et al. 2020).
- `bt.datasets.toy_humann`: a tiny in-memory function table for examples and
  tests.
- `bt.datasets.enzyme`: the ENZYME EC hierarchy (CC BY 4.0), downloaded once
  and cached.
- `bt.datasets.hmp2`: the HMP2 inflammatory bowel disease cohort's pathways,
  species and metadata, one stool metagenome per participant, downloaded once
  and cached.
- A function guide and a tutorial on the HMP2 cohort.
- asv benchmarks for `fn.func_glom`, `io.read_humann` and
  `fn.functional_redundancy`, with their baselines in the docs.
- New runtime dependency: `mudata>=0.4`, for function tables.

### Changed

- `bt.pp.relative` divides each value by its sample total, summed in float64:
  a sample whose total is subnormal no longer gives `inf`, and a float32 table
  is no longer off by about 1e-7.
- `bt.pl.bar` no longer gives a group the grey that marks missing values, so
  with eight or more groups plus missing values the colours change.
- The Coming-from-R page marks the phyloseq calls biotapy does not cover yet
  "not in 0.2".
- The docs build gives each notebook cell up to 5 minutes, enough for the
  function tutorial to download the HMP2 tables on a cold cache.

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
