# Reading and writing data

Readers turn a file format into a `TreeData` that follows the one
[data model](data_model.md) every biotapy function relies on.

No format below records whether its table holds counts or proportions, so
every reader infers `uns["biotapy"]["x_kind"]` from the values: whole numbers
are `"counts"`; otherwise, if every sample with a nonzero total sums to 1
(within `1e-3`), `"relative"`; anything else is `"abundance"`. Functions
that need raw counts check this and refuse proportions.

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
Only the `taxonomy` (or `Taxonomy`) key is read; other observation metadata,
such as `confidence`, is not.

### Sample metadata becomes obs columns

Whatever a table's sample metadata holds becomes columns of `obs`, one column
per key, unchanged except that empty values become `NaN` (BIOM has no null,
so writers store a missing value as an empty string).

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
tree/table alignment every biotapy reader uses. A tree sharing no tip with
the table, or a file that is not valid Newick, raises a `ValueError` naming
`tree=`.

### Writing BIOM

`bt.io.write_biom` writes `X`, taxonomy and sample metadata back out as a BIOM
2.1 HDF5 table by default, or BIOM 1.0 JSON with `fmt="json"`; any other
`fmt` raises a `ValueError`:

```python
bt.io.write_biom(tdata, "table.biom")
```

BIOM has no slot for a tree, layers or embeddings: a TreeData's tree and
everything outside `X`, rank columns and `obs` are not written. Sample
metadata is written as text; missing values are written as empty strings and
read back as `NaN`.

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
`ValueError` naming the argument and the payload it expected, and so does a
file that is not a zip archive at all.

### Matching taxonomy and metadata to the table

Taxonomy and sample metadata are matched to the table by id, with the same
checks for both:

- Ids that are only in the taxonomy or metadata file are ignored: a metadata
  file often covers more samples than one table.
- Table ids missing from the file get `NaN`, and one `UserWarning` gives
  their count.
- A file that shares no id with the table, or that repeats an id, raises a
  `ValueError` naming the argument and example ids.

### Sample metadata

`metadata` reads a QIIME 2 sample-metadata TSV, independent of any artifact:

- The ID column is recognized by its header: `id`, `sampleid`, `sample id`,
  `sample-id`, `featureid`, `feature id` and `feature-id` match
  case-insensitively; the legacy `#SampleID`, `#Sample ID`, `#OTUID`,
  `#OTU ID` and `sample_name` headers match exactly.
- Leading `#`-comment lines and blank rows are skipped.
- A header that repeats a column name raises a `ValueError` naming it.
- An optional `#q2:types` row declares each column `categorical` or
  `numeric`; without it, a column becomes numeric only when every value it
  holds parses as a number, matching QIIME 2's own type inference. A types
  row shorter than the header is padded, and an empty type cell means
  "infer". A column declared `numeric` that holds a non-numeric value (such
  as `1,000` or `thirty`) raises a `ValueError` naming the column and the
  values, as QIIME 2 itself does. Missing values become `NaN`; text columns
  are NaN-backed strings, so a column that is empty for every sample of the
  table still saves to h5ad/h5td.
- The file is read as `utf-8-sig`, so a BOM added by Excel does not break the
  ID header.
- `obs` is aligned to the table's sample order, as described above.

## DADA2

`bt.io.read_dada2` reads a DADA2 sequence table (CSV or TSV) as written by R's
`write.csv`/`write.table`, with optional taxonomy and a tree. Any `.csv`
suffix means CSV, so a compressed `seqtab.csv.gz` works too:

```python
tdata = bt.io.read_dada2("seqtab.csv", taxa="taxa.csv", tree="tree.nwk")
```

Only `seqtab` is required.

### Expected orientation

`seqtab` (DADA2's `seqtab`/`seqtab.nochim`) is already samples (rows) x
sequences (columns), so `read_dada2` reads it as-is, with no transpose.
Its column names must be DNA sequences (`A`, `C`, `G`, `T`, `N`); a table
whose columns are sample ids instead - the transposed orientation - raises a
`ValueError` naming the argument. Sample names are read as text exactly as
written: `001` stays `001`, `1e3` stays `1e3`, and a sample named `NA` is a
name, not a missing value.

### ASV naming

Sequences become the feature ids `ASV1..ASVn`, in column order, with the
original sequence kept in `var["sequence"]`. `taxa` (`assignTaxonomy`/
`addSpecies` output: sequences x `Kingdom..Species`) is matched to `seqtab`'s
sequences and its rank columns are normalized exactly like every other
biotapy reader (lowercase names, `NA` becomes `NaN`). The match follows the
same rules as QIIME 2 taxonomy: sequences only in `taxa` are ignored, a
sequence missing from `taxa` gets `NaN` in every rank column with one
`UserWarning` giving the count, and a `taxa` file sharing no sequence with
`seqtab`, or repeating one, raises a `ValueError`.

### Attaching a tree

Pass a Newick file's path as `tree=` to attach a phylogeny in
`vart["phylo"]`. A DADA2 tree's tips are sequences (e.g. from `AlignSeqs`/
`fasttree` on `seqtab`'s column names), so `read_dada2` relabels them to the
matching `ASV1..ASVn` ids before attaching the tree - the tree's tips and the
table's `var_names` must agree for every other biotapy function to see them
as the same features.

Every tip must be a DNA sequence. A tip named by an ASV id raises a
`ValueError` naming `tree=`: biotapy numbers ASVs by `seqtab`'s column order,
which need not match the numbering behind your tree, so matching by id could
attach branches to the wrong features. Sequence tips that are not in
`seqtab` are pruned, with the usual tree/table warning.

## phyloseq

`bt.io.read_phyloseq` reads a phyloseq object saved from R with `saveRDS()`
or `save()`, directly - no R, no rpy2, and phyloseq itself never has to be
installed:

```python
tdata = bt.io.read_phyloseq("ps.rds")
```

### `.rds` vs `.RData`, and `name=`

An `.rds` file (`saveRDS(ps, "ps.rds")`) holds exactly one object, so
`read_phyloseq` reads it directly. An `.RData`/`.rda` file (`save(ps, ...)`)
can hold several objects; `read_phyloseq` reads the one phyloseq object it
finds automatically. If it holds more than one phyloseq object, pass `name=`
with the R variable name to select one - without it, a `ValueError` lists the
names it found. A file with no phyloseq object at all also raises a
`ValueError`.

### What is read

- `otu_table` becomes `X`, transposed exactly once so samples are always rows
  regardless of whether the R object stored `taxa_are_rows=TRUE` or `FALSE`.
- `tax_table` becomes rank columns in `var`, normalized like every other
  biotapy reader.
- `sample_data` becomes columns of `obs`.
- `phy_tree` becomes the phylogeny in `vart["phylo"]`; ape's internal node
  numbering is replaced with the same collision-free `n0, n1, ...` names
  `read_qiime2`/`read_biom` use for Newick trees.

Any of `tax_table`, `sample_data` or `phy_tree` being absent (`NULL` in R)
gives an empty `var`/`obs` or no tree, never an error.

### What is not read

`refseq` (a `Biostrings` `DNAStringSet` of representative sequences) cannot
be parsed by the underlying `rdata` reader yet. If it is populated,
`read_phyloseq` raises a `ValueError` rather than guessing at or silently
dropping the sequences; the message gives the R fix - export the sequences
separately and re-save the object without the slot:

```r
Biostrings::writeXStringSet(refseq(ps), "refseq.fasta")
ps@refseq <- NULL
saveRDS(ps, "ps.rds")
```
