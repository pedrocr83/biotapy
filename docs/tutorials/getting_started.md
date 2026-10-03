---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Getting started

Install with `pip install biotapy` (Python 3.12 or newer). This page downloads
phyloseq's GlobalPatterns dataset once (435 kB) and draws its first ordination:

```{code-cell} ipython3
import biotapy as bt

tdata = bt.datasets.global_patterns()  # 26 samples x 19,216 OTUs, with taxonomy and a tree
bt.tl.beta(tdata, inplace=True)  # Bray-Curtis distances in obsp["braycurtis"]
bt.tl.pcoa(tdata, inplace=True)  # coordinates in obsm["X_pcoa"]
bt.pl.ordination(tdata, color="SampleType");
```

Each `tl` call with `inplace=True` stores its result in the object; `pl` draws what is stored.
Alpha diversity works the same way:

```{code-cell} ipython3
bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)  # obs["alpha_shannon"]
bt.pl.richness(tdata, "shannon", x="SampleType");
```

Next: the [quick tour](quick_tour.md) of the data model, and the
[phyloseq vignette](phyloseq_analysis.md) redone with biotapy.
