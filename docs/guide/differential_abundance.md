# Differential abundance

Differential abundance (DA) asks which features are more abundant in one group of samples than in
another. Every `bt.da` method takes the same arguments and returns the same table, so their answers
can be put side by side.

## One question, one table

A method takes the `obs` column whose effect you want (`group`), the `covariates` to adjust for,
and the `reference` level that the other level is compared with:

```python
import biotapy as bt

tdata = bt.datasets.toy()
table = bt.da.linda(tdata, "group")  # "B vs A": A, the first category, is the reference
table = bt.da.linda(tdata, "group", reference="B")  # "A vs B": every effect changes sign
```

`group` is a categorical, string or bool column with two levels, or a numeric column, whose effect
is then a slope. Covariates are numeric, or categorical with one indicator per level against
their first. There are no formula strings: R and patsy take the alphabetically first level as the
reference and silently drop samples with a missing value. biotapy takes the reference you give, or
the first category, says which way round the effect is in the `contrast` column, and raises on a
missing value, on a constant column, on a group with more than two levels and on collinear covariates.

The table has one row per feature, in `var_names` order:

| Column | Meaning |
|---|---|
| `effect` | log2 fold change of the other level over the reference, or the slope of a numeric group: per standard deviation in LinDA and MaAsLin 3, per unit in ANCOM-BC2, not offered by ALDEx2. ALDEx2's is the median log2 difference between the groups, MaAsLin 3's the coefficient minus the median coefficient |
| `se` | its standard error; NaN for every feature in ALDEx2, which reports none |
| `pvalue` | the method's p-value |
| `qvalue` | Benjamini-Hochberg adjusted p-value, over the features the method tested |
| `direction` | sign of `effect`: -1, 0 or 1 |
| `method` | the method's name |
| `contrast` | `"<level> vs <reference>"`, or the name of a numeric group |

A feature a method cannot test keeps its row, with NaN `effect`, `pvalue` and `qvalue` and `direction` 0 (a NaN `se` alone means
nothing for ALDEx2). Methods never
filter features: filter once, with `bt.pp.filter_features`, before running any method, so every
method tests the same features. Methods need raw counts in `X` and a read in every sample
(`bt.pp.filter_samples(tdata, min_depth=1)` drops empty ones).

## Choosing methods

biotapy has four methods. They model compositional counts differently, so they call different
features, and no one of them is right on every dataset:

| Method | Model | Runs in | Covariates | Numeric `group` | `seed` |
|---|---|---|---|---|---|
| {doc}`bt.da.linda </methods/linda>` | least squares of log2 centred log-ratios, minus the mode of all coefficients | Python | yes | per standard deviation | no |
| {doc}`bt.da.ancombc2 </methods/ancombc2>` | least squares of log counts corrected for each sample's sampling fraction | Python (scikit-bio) | yes | per unit | no |
| {doc}`bt.da.aldex2 </methods/aldex2>` | Welch t-tests on Monte Carlo draws of centred log-ratios | R | no | no | yes |
| {doc}`bt.da.maaslin3 </methods/maaslin3>` | least squares of log2 relative abundance where present, tested against the median | R | yes | per standard deviation | yes |

Each method's page gives its model, the R defaults biotapy changes, and how closely biotapy
matches the R package on real data. Two things to know before you choose:

- ANCOM-BC2 and MaAsLin 3 leave zeros out, so a feature with no read in one group cannot be
  fitted: its row is NaN. LinDA adds 0.5 to every count when the table has a zero, and ALDEx2
  draws around every count, so both test it.
- `reference` changes more than the sign in ANCOM-BC2 (its bias correction is fitted against the
  reference level) and ALDEx2 (its effect's random resampling follows the label order), and
  slightly in MaAsLin 3 (its p-values move by up to 1e-3 with the same `seed`, the calls do not).
  Choose it on the biology, the control or baseline level, before you look at any result.

(da-methods-in-r)=
## Methods that run in R

ALDEx2 and MaAsLin 3 exist only in R, so biotapy calls them there through
[rpy2](https://rpy2.github.io/). They need R, the R package, and biotapy's `r` extra, which builds
rpy2 against that R (Linux and macOS; it fails to install when no R is found):

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("ALDEx2", "maaslin3"))'
pip install 'biotapy[r]'
```

rpy2 is GPL-2.0-or-later and the R packages have their own licences; biotapy itself does not ship
any of them. Without rpy2 or the R package, the call raises an `ImportError` that names what to
install. Each call converts `X` to a dense table once, because rpy2 has no sparse converter.

Both methods draw random numbers in R. `seed` seeds R for the call, and your R session's own
random state is restored afterwards, so the same seed gives the same table. Each R warning raised
during the call, such as ALDEx2's for fewer than 128 `mc_samples`, is re-emitted as a Python
`UserWarning`; an R error is raised as a `RuntimeError`. Repeated `var_names` or `obs_names`
raise: call `adata.var_names_make_unique()` first.

```python
table = bt.da.aldex2(tdata, "group", seed=0)
table = bt.da.maaslin3(tdata, "group", seed=0)
```

## Where methods agree

`bt.da.consensus` puts the tables of several methods side by side and counts, for each feature,
the methods that call it significant (`qvalue < alpha`, strictly) and whether they agree on its
direction:

```python
results = [bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")]
table = bt.da.consensus(results)  # alpha=0.05, min_methods=2
table[table["consensus"]]
```

A feature is a consensus hit when at least `min_methods` methods call it and all of them give it
the same sign. Methods that call it in opposite directions mark a `conflict`, which is never a
consensus. A method that could not test a feature does not count against it: `n_tested` says how
many methods tested each feature. The tables must come from different methods and compare the same
`contrast`, so run every method with the same `group` and `reference`. The
{doc}`consensus page </methods/consensus>` defines every column.

Methods disagree a lot on real data, and the literature is split on what to do about it: Nearing
et al. (2022) recommend a consensus of several methods, Pelto et al. (2025) one simple method. Either
way, choose the methods before you look at their results; trying methods until one finds what you
hoped for is selective reporting, and a consensus table does not protect against it. Methods that
share a model, such as ANCOM-BC and ANCOM-BC2, also agree more often for that reason alone.

## Plotting the consensus

`bt.pl.consensus` draws the consensus table as a dot matrix: one row per feature that at least one
method calls, one column per method. A filled dot is a call, red for a positive effect and blue for
a negative one; a hollow dot is a feature the method tested without calling it; no dot, a feature
the method could not test. Consensus features have bold labels. Rows are sorted by the number of
methods calling them, then by their mean absolute effect, and `top` keeps the first 30:

```python
ax = bt.pl.consensus(table, top=20)
```

Feature ids are often accession numbers; to label rows with a rank, rename the table's index first,
for example `table.rename(index=tdata.var["genus"])`.
