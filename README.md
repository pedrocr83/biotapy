# biotapy

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/pedrocr83/biotapy/test.yaml?branch=master
[badge-docs]: https://app.readthedocs.org/projects/biotapy/badge/

biotapy is a microbiome analysis toolkit for Python, in the style of R's [mia][]
and [phyloseq][], built on [AnnData][]/[TreeData][]. Samples are always rows,
the count matrix stays sparse, the phylogenetic tree sits alongside the data
as a first-class object, and every computation with an R equivalent is checked
against R on real data.

## Status

**biotapy 0.3 is an early release.** The API can still change between minor versions.

What 0.3 does (full signatures in the [API reference][api]):

- **Readers**: `bt.io.read_biom` (BIOM 1.0/2.1), `bt.io.read_qiime2`
  (`.qza` artifacts, no QIIME 2 install needed), `bt.io.read_dada2`
  (CSV/TSV/`.rds` sequence tables), `bt.io.read_phyloseq`
  (`.rds`/`.RData`, no R needed) and `bt.io.read_metaphlan` (MetaPhlAn 3
  and 4 profiles) all read into one `TreeData`. `bt.io.read_humann` and
  `bt.io.read_picrust2` read function tables into a `MuData` with a community
  and a per-taxon modality; `bt.io.read_picrust2_traits` reads PICRUSt2's
  per-ASV gene copy numbers.
- **Writer**: `bt.io.write_biom` writes a BIOM 2.1 or 1.0 table back out.
- **Datasets**: `bt.datasets.toy` and `bt.datasets.toy_humann` (in-memory,
  for examples and tests), `bt.datasets.global_patterns`,
  `bt.datasets.enterotype` and `bt.datasets.esophagus` (phyloseq's example
  datasets), `bt.datasets.hmp2` (the HMP2 inflammatory bowel disease cohort)
  and `bt.datasets.enzyme` (the ENZYME EC hierarchy), downloaded and cached
  on first use.
- **Preprocessing**: `bt.pp.relative`, `bt.pp.tax_glom`,
  `bt.pp.filter_features`, `bt.pp.filter_samples` and `bt.pp.rarefy`, and
  the compositional transforms `bt.pp.clr` and `bt.pp.philr` (checked against
  vegan and philr).
- **Function**: `bt.fn.load_hierarchy`, `bt.fn.func_glom` and `bt.fn.renorm`
  (checked against HUMAnN's own output), `bt.fn.contributions` and
  `bt.fn.functional_redundancy`.
- **Tools**: `bt.tl.alpha`, `bt.tl.beta`, `bt.tl.unifrac`, `bt.tl.pcoa`,
  `bt.tl.nmds` and `bt.tl.permanova`, each checked against R on real data.
- **Differential abundance**: `bt.da.linda` and `bt.da.ancombc2` in Python,
  `bt.da.aldex2` and `bt.da.maaslin3` through R (the `r` extra below), all
  returning one result table and checked against their R packages, and
  `bt.da.consensus`, which reports where the methods agree.
- **Plots**: `bt.pl.bar`, `bt.pl.richness`, `bt.pl.ordination`,
  `bt.pl.scree`, `bt.pl.heatmap`, `bt.pl.contributions` and
  `bt.pl.consensus`.

Next, in 0.4: multi-omics conventions on MuData, leak-free scikit-learn
transformers, a PyTorch loader and an interface for embedding models. See the
[roadmap][roadmap]; no dates are promised.

## Installation

You need Python 3.12 or newer.

```bash
pip install biotapy
```

or, with [uv][]:

```bash
uv add biotapy
```

or, with [pixi][], in a workspace that has Python 3.12 or newer:

```bash
pixi add "python>=3.12"
pixi add --pypi biotapy
```

The development version installs straight from GitHub:
`pip install git+https://github.com/pedrocr83/biotapy.git`.

ALDEx2 and MaAsLin 3 run in R. They need R with the R packages
(`BiocManager::install(c("ALDEx2", "maaslin3"))`) and the `r` extra, which
builds rpy2 against that R (Linux and macOS):

```bash
pip install 'biotapy[r]'
```

rpy2 is GPL-2.0-or-later and the R packages carry their own licences; biotapy
does not ship any of them.

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

`bt.datasets.global_patterns()`, `bt.datasets.enterotype()` and
`bt.datasets.esophagus()` download phyloseq's example datasets once and cache
them on disk with [pooch][]; later calls reuse the cache with no network
access, as long as the hash still matches. Set `BIOTAPY_DATA_DIR` to change
the cache directory:

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
ones diversity and ordination write) and the
[data-model-slots contract][data-model-contract] for the exact rules.

Every function is pure by default: it returns a new object or its result and
leaves the one you passed in unchanged. Only a `tl` function called with
`inplace=True` writes its result into it (`obs`, `obsp`, `obsm` and
`uns["biotapy"]`).

## Coming from R

Every function's docstring names its R equivalent (`phyloseq::tax_glom`,
`mia::agglomerateByRank`, and so on) in its `Notes` section, and the
[Coming from R][coming-from-r] page lists them all, with the phyloseq accessors
that are plain AnnData code. The [phyloseq analysis vignette][vignette] is
redone with biotapy in the tutorials.

## Example data and licensing

`bt.datasets.global_patterns`, `bt.datasets.enterotype` and
`bt.datasets.esophagus` download example data from
[phyloseq's repository][phyloseq-data] at runtime; biotapy ships none of it.
That data stays licensed to phyloseq's authors under AGPL-3.
`bt.datasets.enzyme` downloads ENZYME, distributed by the SIB Swiss Institute
of Bioinformatics under CC BY 4.0. `bt.datasets.hmp2` downloads the HMP2
tables from the [IBDMDB][ibdmdb], which states no licence for them; cite
Lloyd-Price et al. (2019, *Nature* 569:655-662) when you use them. biotapy
itself is [BSD-3-Clause][license].

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
[ibdmdb]: https://ibdmdb.org/
[license]: https://github.com/pedrocr83/biotapy/blob/master/LICENSE
[rules]: https://github.com/pedrocr83/biotapy/blob/master/rules.md
[knowledge]: https://github.com/pedrocr83/biotapy/tree/master/.knowledge
[golden-playbook]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/playbooks/regenerate-golden-files.md
[roadmap]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/roadmap/index.md
[data-model-contract]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/contracts/data-model-slots.md
[uv]: https://github.com/astral-sh/uv
[pixi]: https://pixi.sh/
[issue tracker]: https://github.com/pedrocr83/biotapy/issues
[tests]: https://github.com/pedrocr83/biotapy/actions/workflows/test.yaml
[documentation]: https://biotapy.readthedocs.io
[changelog]: https://biotapy.readthedocs.io/page/changelog.html
[api]: https://biotapy.readthedocs.io/page/api.html
[guide]: https://biotapy.readthedocs.io/page/guide/index.html
[data-model]: https://biotapy.readthedocs.io/page/guide/data_model.html
[biom-format-1004]: https://github.com/biocore/biom-format/pull/1004
[coming-from-r]: https://biotapy.readthedocs.io/page/coming_from_r.html
[vignette]: https://biotapy.readthedocs.io/page/tutorials/phyloseq_analysis.html
