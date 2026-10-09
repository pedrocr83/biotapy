# Machine learning

`bt.ml` holds scikit-learn transformers, so microbiome preprocessing can sit
inside a [Pipeline](https://scikit-learn.org/stable/modules/compose.html) and
be fitted on the training samples of each cross-validation fold only.

## Which steps leak

A step leaks when what it does to one sample depends on other samples. Run
before the data is split, it lets the test samples shape the training data,
and cross-validated scores come out better than they will be on new samples.

| Step | Learns from the samples? | In a pipeline |
|---|---|---|
| Prevalence filter (`bt.ml.PrevalenceFilter`) | yes: which features are common enough | must be inside |
| Relative abundance (`sklearn.preprocessing.Normalizer(norm="l1")`) | no: each sample on its own | either way |
| CLR (`bt.ml.CLR`) | no: each sample on its own | either way |
| Scaling (`sklearn.preprocessing.StandardScaler`) | yes: each feature's mean and spread | must be inside |

`bt.pp.filter_features` on the whole table before cross-validation is the
leaky version of `bt.ml.PrevalenceFilter`: the same rule, fitted on every
sample at once.

## A leak-free pipeline

```python
import biotapy as bt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_validate
from sklearn.pipeline import make_pipeline

tdata = bt.datasets.toy()
pipeline = make_pipeline(
    bt.ml.PrevalenceFilter(min_prevalence=0.1),
    bt.ml.CLR(pseudocount=0.5),
    LogisticRegression(),
)
scores = cross_validate(pipeline, tdata.X, tdata.obs["group"], cv=3)
```

The transformers take what scikit-learn takes: a samples x features array, a
sparse matrix (`tdata.X`) or a DataFrame, whose column names they keep with
`set_output(transform="pandas")`. They do not take an AnnData.

`PrevalenceFilter` keeps a sparse `X` sparse. `CLR` returns a dense array,
because a log-ratio has no zeros, and gives the values of `bt.pp.clr`; like
`bt.pp.clr` it warns when the pseudocount is larger than the smallest non-zero
value, as when the default 0.5 meets relative abundances.

## Relative abundance

scikit-learn already has it: `Normalizer(norm="l1")` divides each sample by
its total, keeps a sparse matrix sparse and leaves an all-zero sample at zero,
as `bt.pp.relative` does. Before `CLR` it changes nothing but the scale the
pseudocount is on.
