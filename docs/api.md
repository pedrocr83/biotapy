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
    io.read_humann
    io.read_metaphlan
    io.read_phyloseq
    io.read_picrust2
    io.read_picrust2_traits
    io.read_qiime2
    io.write_biom
```

## Datasets

```{eval-rst}
.. module:: biotapy.datasets
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    datasets.enterotype
    datasets.enzyme
    datasets.esophagus
    datasets.global_patterns
    datasets.toy
    datasets.toy_humann
```

## Preprocessing

```{eval-rst}
.. module:: biotapy.pp
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    pp.filter_features
    pp.filter_samples
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
    tl.nmds
    tl.pcoa
    tl.permanova
    tl.unifrac
```

## Plots

```{eval-rst}
.. module:: biotapy.pl
.. currentmodule:: biotapy

.. autosummary::
    :toctree: generated

    pl.bar
    pl.heatmap
    pl.ordination
    pl.richness
    pl.scree
```
