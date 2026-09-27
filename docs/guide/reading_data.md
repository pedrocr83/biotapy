# Reading and writing data

Readers turn a file format into a `TreeData` that follows the one
[data model](data_model.md) every biotapy function relies on.

## BIOM

`bt.io.read_biom` reads both BIOM dialects - JSON 1.0 and HDF5 2.1 - through
[biom-format](http://biom-format.org/), which detects the dialect from the file itself:

```python
import biotapy as bt

tdata = bt.io.read_biom("table.biom")
```

### Samples become rows

BIOM stores its matrix as features by samples. `read_biom` transposes it exactly
once, so `X` comes out samples (rows) by features (columns), like every other
biotapy object.

### Taxonomy becomes rank columns

Observation metadata's `taxonomy` entry - a list of ranks in HDF5, or whatever
string or list a producer wrote in JSON - is parsed into the canonical
`kingdom` .. `species` columns in `var`, the same way for every taxonomy
dialect (Greengenes, RESCRIPT, SILVA). A table with no `taxonomy` metadata
gets a `var` with no rank columns at all, rather than columns full of `NaN`.

### Sample metadata becomes obs columns

Whatever a table's sample metadata holds becomes columns of `obs`, one column
per key, unchanged.

### Ids are always strings

BIOM allows unquoted integer ids in its JSON dialect. `read_biom` casts every
sample and observation id to `str`, so `obs_names`/`var_names` are always
strings regardless of what the file stored.

### Duplicate ids

biom-format validates ids itself: a table with a repeated sample or
observation id raises `biom.exception.TableException` before biotapy ever
sees the table.

### Attaching a tree

Pass a Newick file's path as `tree=` to attach a phylogeny in `vart["phylo"]`.
If the tree's tips and the table's features disagree, only the shared
features are kept and one `UserWarning` names both counts - the same
tree/table alignment every biotapy reader uses.

### Writing BIOM

`bt.io.write_biom` writes `X`, taxonomy and sample metadata back out as a BIOM
2.1 HDF5 table by default, or BIOM 1.0 JSON with `fmt="json"`:

```python
bt.io.write_biom(tdata, "table.biom")
```

BIOM has no slot for a tree, layers or embeddings: a TreeData's tree and
everything outside `X`, rank columns and `obs` are not written. Sample
metadata is written as text; missing values become empty strings and read
back as `""`, not NaN.

Rank columns are written as prefixed values (`k__`, `p__`, ..., `g__`) rather
than bare strings, even where a rank is missing. BIOM's HDF5 reader drops
empty taxonomy entries on read, which would otherwise shift every rank after
a missing one out of place; a bare prefix like `g__` stays truthy, so
`read_biom` can map it back to the right rank by its letter.

## QIIME 2

`bt.io.read_qiime2` reads QIIME 2 artifacts (`.qza`) directly, without a
QIIME 2 install: a `.qza` is just a zip file, and biotapy extracts the
payload it needs from it.

```python
tdata = bt.io.read_qiime2(
    "table.qza", taxonomy="taxonomy.qza", tree="tree.qza", metadata="sample-metadata.tsv"
)
```

Only `table` is required. Each argument reads one artifact:

- `table`: a `FeatureTable[Frequency]` artifact (`data/feature-table.biom`,
  BIOM 2.1 HDF5), read the same way as `read_biom`.
- `taxonomy`: a `FeatureData[Taxonomy]` artifact (`data/taxonomy.tsv`); its
  `Taxon` column is parsed into rank columns exactly like BIOM taxonomy, and
  an optional `Confidence` column becomes a `confidence` column in `var`.
  Passing `taxonomy` replaces any taxonomy already in the table.
- `tree`: a `Phylogeny[Rooted]` or `Phylogeny[Unrooted]` artifact
  (`data/tree.nwk`), attached the same way as `read_biom`'s `tree` argument.

An artifact is recognized by the payload file it holds, not by its
`metadata.yaml`, so any `.qza` holding the expected payload works. Passing
the wrong artifact - for example a taxonomy artifact as `table` - raises a
`ValueError` naming the argument and the payload it expected.

### Sample metadata

`metadata` reads a QIIME 2 sample-metadata TSV, independent of any artifact:

- The ID column is recognized by its header: `id`, `sampleid`, `sample id`,
  `sample-id`, `featureid`, `feature id` and `feature-id` match
  case-insensitively; the legacy `#SampleID`, `#Sample ID`, `#OTUID`,
  `#OTU ID` and `sample_name` headers match exactly.
- Leading `#`-comment lines and blank rows are skipped.
- An optional `#q2:types` row declares each column `categorical` or
  `numeric`; without it, a column becomes numeric only when every value it
  holds parses as a number, matching QIIME 2's own type inference. A types
  row shorter than the header is padded, and an empty type cell means
  "infer". Missing values become `NaN`.
- The file is read as `utf-8-sig`, so a BOM added by Excel does not break the
  ID header.
- Samples in the metadata that are not in the table are ignored; `obs` is
  aligned to the table's sample order.

## DADA2

`bt.io.read_dada2` reads a DADA2 sequence table (CSV or TSV) as written by R's
`write.csv`/`write.table`, with optional taxonomy and a tree:

```python
tdata = bt.io.read_dada2("seqtab.csv", "taxa.csv", tree="tree.nwk")
```

Only `seqtab` is required.

### Expected orientation

`seqtab` (DADA2's `seqtab`/`seqtab.nochim`) is already samples (rows) x
sequences (columns), so `read_dada2` reads it as-is, with no transpose.
Its column names must be DNA sequences (`A`, `C`, `G`, `T`, `N`); a table
whose columns are sample ids instead - the transposed orientation - raises a
`ValueError` naming the argument.

### ASV naming

Sequences become the feature ids `ASV1..ASVn`, in column order, with the
original sequence kept in `var["sequence"]`. `taxa` (`assignTaxonomy`/
`addSpecies` output: sequences x `Kingdom..Species`) is matched to `seqtab`'s
sequences and its rank columns are normalized exactly like every other
biotapy reader (lowercase names, `NA` becomes `NaN`); a sequence missing from
`taxa` gets `NaN` in every rank column.

### Attaching a tree

Pass a Newick file's path as `tree=` to attach a phylogeny in
`vart["phylo"]`. A DADA2 tree's tips are sequences (e.g. from `AlignSeqs`/
`fasttree` on `seqtab`'s column names), so `read_dada2` relabels them to the
matching `ASV1..ASVn` ids before attaching the tree - the tree's tips and the
table's `var_names` must agree for every other biotapy function to see them
as the same features.
