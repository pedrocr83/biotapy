# API

Public functions are listed here as they ship, from Phase 1 onward.

## Input and output

```{eval-rst}
.. module:: biotapy.io
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    io.read_biom
    io.read_dada2
    io.read_h5mu
    io.read_humann
    io.read_metaphlan
    io.read_phyloseq
    io.read_picrust2
    io.read_picrust2_traits
    io.read_qiime2
    io.to_mudata
    io.write_biom
    io.write_h5mu
```

## Datasets

```{eval-rst}
.. module:: biotapy.datasets
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    datasets.biocrust
    datasets.enterotype
    datasets.enzyme
    datasets.esophagus
    datasets.global_patterns
    datasets.hmp2
    datasets.toy
    datasets.toy_humann
```

## Preprocessing

```{eval-rst}
.. module:: biotapy.pp
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    pp.clr
    pp.filter_features
    pp.filter_samples
    pp.philr
    pp.rarefy
    pp.relative
    pp.tax_glom
```

## Function

```{eval-rst}
.. module:: biotapy.fn
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    fn.contributions
    fn.func_glom
    fn.functional_redundancy
    fn.load_hierarchy
    fn.renorm
```

## Tools

```{eval-rst}
.. module:: biotapy.tl
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    tl.alpha
    tl.beta
    tl.mmvec
    tl.nmds
    tl.pcoa
    tl.permanova
    tl.unifrac
```

## Differential abundance

```{eval-rst}
.. module:: biotapy.da
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    da.aldex2
    da.ancombc2
    da.consensus
    da.linda
    da.maaslin3
```

## Machine learning

```{eval-rst}
.. module:: biotapy.ml
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    ml.CLR
    ml.PrevalenceFilter
    ml.embed
    ml.to_torch
```

## Plots

```{eval-rst}
.. module:: biotapy.pl
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    pl.bar
    pl.consensus
    pl.contributions
    pl.heatmap
    pl.ordination
    pl.richness
    pl.scree
```
