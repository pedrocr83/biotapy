# biotapy

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/pedrocr83/biotapy/test.yaml?branch=master
[badge-docs]: https://app.readthedocs.org/projects/biotapy/badge/

biotapy is a microbiome analysis toolkit for Python, in the style of R's [mia][]
and [phyloseq][], built on [AnnData][]/[TreeData][]. Samples are always rows,
the count matrix stays sparse, the phylogenetic tree sits alongside the data
as a first-class object, and every function with an R equivalent is checked
against R on real data.

## Status

**biotapy is pre-release.** The API can still change without notice.

What works today (full signatures in the [API reference][api]):

- **Readers**: `bt.io.read_biom` (BIOM 1.0/2.1), `bt.io.read_qiime2`
  (`.qza` artifacts, no QIIME 2 install needed), `bt.io.read_dada2`
  (CSV/TSV/`.rds` sequence tables) and `bt.io.read_phyloseq`
  (`.rds`/`.RData`, no R needed) all read into one `TreeData`.
- **Writer**: `bt.io.write_biom` writes a BIOM 2.1 or 1.0 table back out.
- **Datasets**: `bt.datasets.toy` (in-memory, for examples and tests),
  `bt.datasets.global_patterns` and `bt.datasets.enterotype`
  (phyloseq's example datasets, downloaded and cached on first use).
- **Preprocessing**: `bt.pp.relative` (per-sample relative abundance)
  and `bt.pp.tax_glom` (aggregate to a taxonomic rank).

What comes next, in release 0.1: filtering (`filter_features`/`filter_samples`),
rarefaction, alpha and beta diversity, UniFrac, PCoA/NMDS ordination,
PERMANOVA, plots, and a generated Coming-from-R table. See the
[Phase 1 roadmap][roadmap] for the full list; no dates are promised.

## Installation

You need Python 3.12 or newer.

biotapy is not functional on PyPI yet: `biotapy` 0.0.1 there is a name
placeholder with no public functions, reserved ahead of the first real
release. The first functional PyPI release will be 0.1. Until then, install
the development version straight from GitHub:

```bash
pip install git+https://github.com/pedrocr83/biotapy.git
```

or, with [uv][]:

```bash
uv add git+https://github.com/pedrocr83/biotapy.git
```

<!--
Once 0.1 is on PyPI:

```bash
uv add biotapy
```

```bash
pip install biotapy
```
-->

On Python 3.14, the `biom-format` dependency has no wheels yet, so it is built
from source and needs a C compiler until biom-format publishes 3.14 wheels
([biocore/biom-format#1004][biom-format-1004]).

## Quick start

Load the bundled toy dataset, add relative abundance, and aggregate to phylum:

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.shape  # (6, 8): 6 samples, 8 features

rel = bt.pp.relative(tdata)
rel.layers["relative"][0].sum()  # 1.0, X itself is untouched

glom = bt.pp.tax_glom(tdata, "phylum")
glom.n_vars  # 3: features sharing a lineage down to phylum are merged
```

### Reading your own data

Each reader returns a `TreeData` following the same [data model](#the-data-model).
Only the first argument is required; the keyword arguments attach taxonomy,
sample metadata and a phylogenetic tree where the format supports them.

```python
# BIOM 1.0 (JSON) or 2.1 (HDF5), dialect detected automatically
tdata = bt.io.read_biom("table.biom", tree="tree.nwk")

# QIIME 2 artifacts (.qza), no QIIME 2 install required
tdata = bt.io.read_qiime2(
    "table.qza", taxonomy="taxonomy.qza", tree="tree.qza", metadata="sample-metadata.tsv"
)

# DADA2 sequence tables: CSV, TSV or R's .rds, with optional taxonomy and tree
tdata = bt.io.read_dada2("seqtab.rds", taxa="taxa.rds", tree="tree.nwk")

# phyloseq objects saved from R with saveRDS()/save() - no R, no rpy2
tdata = bt.io.read_phyloseq("ps.rds")
# an .RData file can hold several phyloseq objects; pick one by its R variable name
tdata = bt.io.read_phyloseq("ps.RData", name="ps_16S")
```

Write back out to BIOM with `bt.io.write_biom(tdata, "table.biom")`.

### Example datasets

`bt.datasets.global_patterns()` and `bt.datasets.enterotype()` download
phyloseq's example datasets once and cache them on disk with [pooch][]; later
calls reuse the cache with no network access, as long as the hash still
matches. Set `BIOTAPY_DATA_DIR` to change the cache directory:

```python
import biotapy as bt

gp = bt.datasets.global_patterns()  # 26 samples, 19,216 OTUs
gp.shape
```

### Saving your work

A `TreeData` round-trips through [treedata][]'s own format, tree included:

```python
tdata.write_h5td("tdata.h5td")

import treedata as td

back = td.read_h5td("tdata.h5td")
```

Dropping down to plain AnnData (no tree) works the usual way too:
`tdata.write_h5ad("tdata.h5ad")` / `anndata.read_h5ad("tdata.h5ad")`.

## The data model

Every biotapy function reads and writes the same object, a `TreeData` (an
AnnData that also carries a phylogeny):

| Slot | Holds |
| --- | --- |
| `X` | Samples x features, sparse CSR |
| `obs` | Sample metadata |
| `var` | Taxonomy, one column per rank: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` |
| `vart["phylo"]` | The phylogenetic tree, leaves named after your features |
| `layers["relative"]` | Relative abundance, added by `pp.relative` |
| `uns["biotapy"]` | `x_kind` (what `X` currently holds: counts, relative, ...) and `provenance` (which biotapy functions produced this object, in order) |

See the [data model guide][data-model] for the full slot list (including the
ones diversity and ordination will use once they ship) and the
[data-model-slots contract][data-model-contract] for the exact rules.

Every function is **pure by default**: it returns a new object and never
changes the one you passed in.

## Coming from R

Every function's docstring names its R equivalent (`phyloseq::tax_glom`,
`mia::agglomerateByRank`, and so on) in its `Notes` section. A generated,
searchable Coming-from-R table arrives with release 0.1.

## Example data and licensing

`bt.datasets.global_patterns` and `bt.datasets.enterotype` download example
data from [phyloseq's repository][phyloseq-data] at runtime; biotapy ships
none of it. That data stays licensed to phyloseq's authors under AGPL-3.
biotapy itself is [BSD-3-Clause][license].

## Documentation

Full documentation, including the [user guide][guide] and [API reference][api],
is at [biotapy.readthedocs.io][documentation].

## Development

```bash
uv sync --group dev --group test --group doc
uvx prek run --all-files
uv run --group test pytest
```

Network and golden tests (downloads data, compares against R) are excluded by
default; run them explicitly with:

```bash
uv run --group test pytest -m "network or golden"
```

Regenerating the R golden files themselves needs the pinned R container - see
the [regenerate-golden-files playbook][golden-playbook].

All contributions follow [rules.md][rules] and the
[knowledge bundle][knowledge] (decisions, contracts, roadmap).

## Release notes

See the [changelog][].

## Contact

Questions, bug reports and feature requests all go to the [issue tracker][].

## Citation

> t.b.a

[mia]: https://microbiome.github.io/mia/
[phyloseq]: https://joey711.github.io/phyloseq/
[anndata]: https://anndata.readthedocs.io/
[treedata]: https://treedata.readthedocs.io/
[pooch]: https://www.fatiando.org/pooch/
[phyloseq-data]: https://github.com/joey711/phyloseq/tree/master/data
[license]: https://github.com/pedrocr83/biotapy/blob/master/LICENSE
[rules]: https://github.com/pedrocr83/biotapy/blob/master/rules.md
[knowledge]: https://github.com/pedrocr83/biotapy/tree/master/.knowledge
[golden-playbook]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/playbooks/regenerate-golden-files.md
[roadmap]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/roadmap/phase-1-core.md
[data-model-contract]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/contracts/data-model-slots.md
[uv]: https://github.com/astral-sh/uv
[issue tracker]: https://github.com/pedrocr83/biotapy/issues
[tests]: https://github.com/pedrocr83/biotapy/actions/workflows/test.yaml
[documentation]: https://biotapy.readthedocs.io
[changelog]: https://biotapy.readthedocs.io/page/changelog.html
[api]: https://biotapy.readthedocs.io/page/api.html
[guide]: https://biotapy.readthedocs.io/page/guide/index.html
[data-model]: https://biotapy.readthedocs.io/page/guide/data_model.html
[biom-format-1004]: https://github.com/biocore/biom-format/pull/1004
