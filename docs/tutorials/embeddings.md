# Embedding samples with MGM

This tutorial runs one pretrained microbiome model end to end: GlobalPatterns' genus profiles go
through MGM, each sample comes back as a vector of 256 numbers, and a nearest-neighbour check
shows that the vectors keep samples from the same environment together.

:::{note}
This page is not run when the documentation is built. MGM needs PyTorch and transformers (the
extra `mgm`), which the documentation build does not install. The code below is run on every
pull request instead, in the `ml-extras` job of biotapy's continuous integration, which checks
that it gives the outputs quoted here (`tests/ml/test_mgm.py`).
:::

```bash
pip install 'biotapy[mgm]'
```

## Embedding the samples

MGM reads genera, so the operational taxonomic units are merged to genus first. The first call
downloads MGM's pretrained model once (33 MB) into the cache `bt.datasets` uses.

```python
import biotapy as bt

genera = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
bt.ml.embed(genera, "mgm", inplace=True)
embedding = genera.obsm["X_mgm"]
```

```text
UserWarning: mgm leaves out 96 of 996 features: 0 without a genus and 96 whose genus is not one of MGM's (4-29, 4041AA30, A17, Aquamonas, Arctic95A-2, ...)
```

MGM's vocabulary holds 9,665 genera; 96 of GlobalPatterns' 996 genus-level features name a genus
outside it and are left out. Most are environmental lineages known only by a clone name (`4-29`,
`BD2-13`) or *Candidatus* genera, which GlobalPatterns' Greengenes taxonomy writes as one word
(`CandidatusPelagibacter`) where MGM's vocabulary has `Candidatus_Pelagibacter`; biotapy reads
each name as MGM's own code does. Every sample keeps genera MGM knows, so none is embedded from
an empty profile.

```python
embedding.shape, embedding.dtype
```

```text
((26, 256), dtype('float32'))
```

## What the embedding keeps

MGM never saw GlobalPatterns' sample types. For each sample, the nearest other sample in the
embedding, by cosine distance, should still come from the same environment:

```python
from sklearn.neighbors import NearestNeighbors

nearest = NearestNeighbors(n_neighbors=2, metric="cosine").fit(embedding).kneighbors(embedding)[1][:, 1]
sample_type = genera.obs["SampleType"].to_numpy()
same_type = int((sample_type[nearest] == sample_type).sum())
same_type
```

```text
26
```

For all 26 samples, the nearest neighbour is a sample of the same type: feces next to feces, soil
next to soil, the three mock communities together. The embedding is a starting point for any
model that takes a fixed-length vector per sample; `genera.obsm["X_mgm"]` travels with the
table.

## Time and memory

In one unpinned laptop run, not checked by CI, on a CPU with 8 threads (torch 2.14.1+cpu,
transformers 5.19.0), the first call took 7 s, most of it importing torch and loading the model,
and a second call 1 s for the 26 samples: at least 26 samples a second, in line with the guide's
"a few tens of samples a second". MGM runs one sample at a time and needs no memory beyond the
model's.

## More

The {doc}`machine learning guide </guide/machine_learning>` says how MGM reads a table's genera,
which features it leaves out, how the embedding is pooled, how to cite MGM and its licence, and
how a package adds its own model as a plugin.
