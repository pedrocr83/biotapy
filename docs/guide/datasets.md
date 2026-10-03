# Example datasets

`biotapy.datasets` ships five example datasets. Four return the [data model](data_model.md)
every biotapy function relies on; `toy_humann` returns a `MuData` function
table.

## `toy`

`bt.datasets.toy()` is a tiny 6-sample x 8-feature dataset built in memory -
no download, no network. It carries taxonomy (`kingdom` through `genus`,
with one feature unassigned at `genus`) and a phylogeny, and is what every
docstring `Examples` section and most tests use:

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.shape  # (6, 8)
```

## `toy_humann`

`bt.datasets.toy_humann()` is the toy samples' gene families as HUMAnN writes
them after regrouping to EC numbers, built in memory: a `MuData` whose
`"function"` modality holds `UNMAPPED`, `UNGROUPED` and four EC numbers in
RPK, and whose `"function_by_taxon"` modality holds their seven strata.
`bt.fn` docstrings use it:

```python
import biotapy as bt

mdata = bt.datasets.toy_humann()
mdata["function"].shape  # (6, 6)
```

## `global_patterns`, `enterotype` and `esophagus`

`bt.datasets.global_patterns()`, `bt.datasets.enterotype()` and
`bt.datasets.esophagus()` are three well-known phyloseq example datasets,
read through `bt.io.read_phyloseq`:

- **GlobalPatterns**: 26 samples from 9 environments, 19,216 OTUs, with
  taxonomy and a tree; counts in `X`.
- **enterotype**: 280 gut samples, 553 genera, no tree; relative abundances
  in `X` (`uns["biotapy"]["x_kind"] == "relative"`).
- **esophagus**: 3 esophageal biopsies (samples `B`, `C`, `D`), 58 OTUs, with a tree but no
  taxonomy or sample data; counts in `X`. biotapy's UniFrac golden tests use it.

```python
import biotapy as bt

tdata = bt.datasets.global_patterns()
```

## Caching with pooch

Each of `global_patterns`, `enterotype` and `esophagus` is downloaded once from
[phyloseq's repository][phyloseq-data] - a single pinned commit, so the
`.RData` files' SHA-256 hashes stay valid - and cached on disk with
[pooch](https://www.fatiando.org/pooch/). A later call re-hashes the cached
file and, as long as the hash still matches, reads it straight from disk
with no network access at all.

Set `BIOTAPY_DATA_DIR` to change the cache directory; the default is a
per-user cache directory (`pooch.os_cache("biotapy")`):

```bash
BIOTAPY_DATA_DIR=/path/to/cache python my_script.py
```

## Licensing

`global_patterns`, `enterotype` and `esophagus` download data from phyloseq's repository
at runtime; biotapy ships none of it. That data stays licensed to
phyloseq's authors under AGPL-3. biotapy itself is
[BSD-3-Clause](https://github.com/pedrocr83/biotapy/blob/master/LICENSE).

[phyloseq-data]: https://github.com/joey711/phyloseq/tree/master/data
