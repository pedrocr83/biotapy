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

# Leak-free cross-validation

Cross-validation estimates how well a model does on samples it has not seen. A preprocessing
step that learns from the samples, such as choosing which features to keep, leaks when it is
fitted on every sample before the split: the test samples then shape what the model is trained
on, and the score comes out better than the model will do on new samples. This tutorial measures
that on a real cohort, with `bt.ml`'s transformers inside a scikit-learn pipeline, so that every
step is fitted on the training samples of each fold only.

:::{note}
The HMP2 tables are downloaded from the [IBDMDB](https://ibdmdb.org/) on first use (23 MB) and
cached. The IBDMDB states no licence for them; biotapy ships none of them. Cite the study when
you use them: Lloyd-Price J et al. (2019) Multi-omics of the gut microbial ecosystem in
inflammatory bowel diseases. *Nature* 569:655-662.
:::

```{code-cell} ipython3
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline

import biotapy as bt
```

## The question

Does the species profile of a participant's first stool sample tell inflammatory bowel disease
(Crohn's disease or ulcerative colitis) from non-IBD? `bt.datasets.hmp2()` holds the HMP2
cohort's MetaPhlAn 3 species for 130 participants.

```{code-cell} ipython3
taxa = bt.datasets.hmp2()["taxa"]
taxa.obs["diagnosis"].value_counts()
```

```{code-cell} ipython3
ibd = (taxa.obs["diagnosis"] != "nonIBD").to_numpy()
```

MetaPhlAn's values are relative abundances: each sample sums to 1. `bt.ml.CLR`'s default
pseudocount of 0.5 would swamp them, so it is set on their scale, below the smallest non-zero
value:

```{code-cell} ipython3
f"{taxa.X.data.min():.2g}"
```

## The pipeline

The prevalence filter keeps the species present in at least 10% of the samples, the centred
log-ratio puts them on a log scale, and a logistic regression predicts IBD. Five-fold
cross-validation, stratified so every fold keeps the share of non-IBD participants, is repeated
five times with different splits, and the score is the mean ROC AUC over the 25 test folds:
0.5 is a coin toss, 1 a perfect ranking.

```{code-cell} ipython3
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)


def mean_auc(model, X, y):
    return float(cross_val_score(model, X, y, cv=cv, scoring="roc_auc").mean())


pipeline = make_pipeline(
    bt.ml.PrevalenceFilter(min_prevalence=0.1),
    bt.ml.CLR(pseudocount=1e-6),
    LogisticRegression(max_iter=5000),
)
inside = mean_auc(pipeline, taxa.X, ibd)
round(inside, 3)
```

Close to a coin toss: with this model, a first stool sample's species hardly separate IBD from
non-IBD participants in this cohort. That is the honest answer, and the one the leaky versions
below would have hidden.

## The prevalence filter outside the pipeline

The common mistake is to filter the whole table first, with `bt.pp.filter_features`, and
cross-validate only the model:

```{code-cell} ipython3
filtered = bt.pp.filter_features(taxa, min_prevalence=0.1)
outside = mean_auc(make_pipeline(bt.ml.CLR(pseudocount=1e-6), LogisticRegression(max_iter=5000)), filtered.X, ibd)
round(outside, 3)
```

The prevalence filter fitted inside the pipeline scores 0.527 and the one fitted on every
sample 0.526: no measurable leak. The filter never looks at the labels, so here, fitting it on the
test samples moves the score by 0.001; it only changes which rare species the model gets. A step
fitted on all samples without the labels can still bias cross-validation in general, so keep it
inside the pipeline all the same: it costs nothing, and then no step is fitted on a test sample
without anyone having to argue that this one is harmless.

## A step that reads the labels

Choosing features by how well they separate the groups is different. `SelectKBest` keeps the 20
species whose CLR values differ most between IBD and non-IBD (an F-test). Inside the pipeline,
it chooses them from the training samples of each fold; outside, from all 130 samples, test
samples included:

```{code-cell} ipython3
selecting = make_pipeline(
    bt.ml.PrevalenceFilter(min_prevalence=0.1),
    bt.ml.CLR(pseudocount=1e-6),
    SelectKBest(f_classif, k=20),
    LogisticRegression(max_iter=5000),
)
select_inside = mean_auc(selecting, taxa.X, ibd)
clr = bt.ml.CLR(pseudocount=1e-6).fit_transform(filtered.X)
select_outside = mean_auc(LogisticRegression(max_iter=5000), SelectKBest(f_classif, k=20).fit_transform(clr, ibd), ibd)
round(select_inside, 3), round(select_outside, 3)
```

Selected inside the pipeline, the score is 0.558; selected on every sample, 0.705. The second
number describes a classifier that does not exist: its 20 species were chosen because they
separate the very samples it is then tested on.

## Labels that mean nothing

The plainest check is to shuffle the labels, so that no species can predict them, and run both
versions again. Ten shuffles, each cross-validated as above:

```{code-cell} ipython3
rng = np.random.default_rng(0)
shuffles = [rng.permutation(ibd) for _ in range(10)]
shuffled_scores = [mean_auc(selecting, taxa.X, labels) for labels in shuffles]
shuffled_inside = float(np.mean(shuffled_scores))
shuffled_outside = float(
    np.mean(
        [
            mean_auc(LogisticRegression(max_iter=5000), SelectKBest(f_classif, k=20).fit_transform(clr, labels), labels)
            for labels in shuffles
        ]
    )
)
round(shuffled_inside, 3), round(shuffled_outside, 3), round(min(shuffled_scores), 3), round(max(shuffled_scores), 3)
```

On shuffled labels the pipeline scores 0.544 on average, near a coin toss, as it should. The
version that selects species on every sample scores 0.767: higher than it scored on the real
labels. The leak manufactures a signal from noise.

The ten shuffled-label scores of the pipeline range from 0.484 to 0.669, and the honest scores of
0.527 and 0.558 fall inside that range: on this cohort they cannot be told apart from scores on
shuffled labels.

```{code-cell} ipython3
pd.DataFrame(
    {"inside the pipeline": [inside, select_inside, shuffled_inside], "on every sample": [outside, select_outside, shuffled_outside]},
    index=["prevalence filter", "SelectKBest, real labels", "SelectKBest, shuffled labels"],
).round(3)
```

## What goes inside the pipeline

Every step that learns from more than one sample: which features to keep, by prevalence or by
their relation to the labels, and how to scale them. Steps that transform each sample on its own,
such as relative abundance and CLR, give the same values either way; putting them inside keeps the
rule simple. The [machine learning guide](../guide/machine_learning.md#which-steps-leak) lists
which `bt.ml` and scikit-learn steps learn from the samples.
