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

# Microbes and metabolites in a desert biocrust

This tutorial puts two data types measured on the same samples into one `MuData`, then asks
which metabolites go with which microbe, with mmvec (Morton et al. 2019). The data are the
example from mmvec's own repository: a desert biological soil crust from four successional stages
(`early`, `earlymid`, `latemid`, `late`), sampled at five times after wetting, with its microbes
counted and its metabolites measured (Swenson et al. 2018).

:::{note}
The two tables are downloaded from [mmvec's repository](https://github.com/biocore/mmvec/tree/master/examples/soils)
on first use (135 KB) and cached; they are distributed there under its BSD-3-Clause licence.
Cite the mmvec paper when you use them: Morton JT et al. (2019) Learning representations of
microbe-metabolite interactions. *Nature Methods* 16:1306-1314. The data come from Swenson TL et
al. (2018) Linking soil biology and chemistry in biological soil crust using isolate
exometabolomics. *Nature Communications* 9:19.
:::

```{code-cell} ipython3
import biotapy as bt
```

## Two tables, one MuData

The microbe table has 20 samples and the metabolite table 19; `bt.datasets.biocrust()` keeps the
19 both have, as `bt.io.to_mudata` does for your own tables. Each data type is a modality, named
as biotapy's {doc}`multi-omics guide </guide/multiomics>` names them.

```{code-cell} ipython3
mdata = bt.datasets.biocrust()
mdata
```

## Which metabolites go with which microbe

`bt.tl.mmvec` learns, from the samples, how likely each metabolite is given each microbe. Its
result has one row per microbe and one column per metabolite, holding log probabilities centred
so each row sums to 0: within a row, a large value is a metabolite that tends to be abundant
where that microbe is. The fit starts from random values, so `seed` fixes it.

```{code-cell} ipython3
ranks = bt.tl.mmvec(mdata, seed=0)
ranks.shape
```

The most abundant microbe is `rplo 1 (Cyanobacteria)`, which mmvec's example treats as the
cyanobacterium *Microcoleus vaginatus*. Its ten highest-ranked metabolites:

```{code-cell} ipython3
cyanobacterium = ranks.loc["rplo 1 (Cyanobacteria)"].sort_values(ascending=False)
cyanobacterium.head(10)
```

## A check against mmvec's own example

mmvec's repository checks its fit of these data against 13 metabolites it lists for this
microbe: every one should rank above zero in its row. biotapy fits mmvec through scikit-bio
rather than mmvec's own TensorFlow code, so the same check says whether the answer holds:

```{code-cell} ipython3
microcoleus = {
    "(3-methyladenine)", "7-methyladenine", "4-guanidinobutanoate", "uracil", "xanthine", "hypoxanthine",
    "(N6-acetyl-lysine)", "cytosine", "N-acetylornithine", "succinate", "adenosine", "guanine", "adenine",
}
above_zero = int((cyanobacterium[list(microcoleus)] > 0).sum())
above_zero
```

All 13 of the metabolites mmvec's example lists rank above zero for the cyanobacterium, as in
mmvec's own fit, although biotapy's fit differs in its optimiser and in its number of
dimensions (scikit-bio's default of 3, where the example used 1).

## Reading the result

- A rank compares metabolites within one microbe's row. It is not a correlation, and a low rank
  means no association rather than a negative one.
- 19 samples is few: mmvec's README says studies this small need careful tuning of the number
  of dimensions and the priors, which scikit-bio's `mmvec` exposes and `bt.tl.mmvec` leaves at
  their defaults.
- The {doc}`multi-omics guide </guide/multiomics>` says how to build the `MuData` from your own
  tables and what `bt.tl.mmvec` refuses.
