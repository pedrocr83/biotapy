# Machine learning

`bt.ml` holds scikit-learn transformers, so microbiome preprocessing can sit
inside a [Pipeline](https://scikit-learn.org/stable/modules/compose.html) and
be fitted on the training samples of each cross-validation fold only, turns
a table into a PyTorch dataset, and embeds samples with pretrained models.

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
sample at once. The {doc}`leak-free cross-validation tutorial
</tutorials/leak_free_cv>` measures both versions on the HMP2 cohort: there
the prevalence filter leaks almost nothing, while a step that chooses features
by the labels inflates the score even on shuffled labels.

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

## PyTorch

`bt.ml.to_torch` turns a table into a PyTorch
[dataset](https://docs.pytorch.org/docs/stable/data.html) with one item per
sample, for a `DataLoader` to batch and shuffle. It needs the extra `torch`:

```bash
pip install 'biotapy[torch]'
```

```python
import biotapy as bt
from torch.utils.data import DataLoader

tdata = bt.pp.clr(bt.datasets.toy())
dataset = bt.ml.to_torch(tdata, label_key="group", layer="clr")
for features, labels in DataLoader(dataset, batch_size=32, shuffle=True):
    ...  # features: float32, batch x features; labels: int64 codes
```

An item is a sample's features as a float32 tensor, or with `label_key` the
pair `(features, label)`. A category, string or bool column gives int64 codes
in category order, for a classifier; a numeric column gives float32, for a
regression. A missing label raises: drop those samples first.

Rows are densified one at a time, as they are read, so a sparse table is never
dense in full. The dataset holds the table without copying it and reads the
labels once, so do not modify the AnnData while you use the dataset.

`to_torch` neither splits nor fits anything. Build one dataset over the whole
table and split it with `torch.utils.data.Subset`: label codes are per
dataset, and anndata drops a category a subset lacks, so a dataset built per
split recodes the classes when a split misses one. A step that learns from the
samples, like the prevalence filter, is fitted on the training samples only and
applied to both splits:

```python
from torch.utils.data import Subset

train, test = [0, 1, 3, 4], [2, 5]
keep = bt.ml.PrevalenceFilter(min_prevalence=0.1).fit(tdata[train].X).get_support()
dataset = bt.ml.to_torch(tdata[:, keep], label_key="group")
train_set, test_set = Subset(dataset, train), Subset(dataset, test)
```

torch publishes no wheel for Intel macOS, so the extra does not install there.

On Linux, pip installs PyPI's torch, which brings CUDA libraries. For a
CPU-only torch, install it from PyTorch's CPU index first:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install 'biotapy[torch]'
```

## Embeddings

`bt.ml.embed` runs a pretrained model over every sample and returns one
vector per sample, in `obs` order, or with `inplace=True` stores it in
`obsm["X_<model>"]`:

```python
import biotapy as bt

genera = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
bt.ml.embed(genera, "mgm", inplace=True)
genera.obsm["X_mgm"].shape  # (26, 256)
```

The {doc}`embedding tutorial </tutorials/embeddings>` runs this on
GlobalPatterns and checks what the embedding keeps.

### MGM

`"mgm"` is MGM, the Microbial General Model of Zhang et al. (2026): a GPT-2
with 8 layers and 256 dimensions, pretrained on genus profiles from MGnify,
released under the MIT licence. biotapy ships it as a plugin behind the extra
`mgm`, which installs torch and transformers:

```bash
pip install 'biotapy[mgm]'
```

The first call downloads the pretrained model once: 33 MB, the
`microformer-mgm` 0.5.8 wheel from PyPI, checked against its SHA-256, into the
cache `bt.datasets` uses (`BIOTAPY_DATA_DIR` if set, otherwise pooch's
per-user cache directory). `microformer-mgm` itself is never installed.

MGM reads genera. Each feature's `var["genus"]` is read as MGM reads a
`g__<genus>` column: up to the first character that is not a letter, digit
or underscore, so `Escherichia-Shigella` is MGM's `Escherichia`. Features of
the same genus are summed, so a table at any level works, and
`bt.pp.tax_glom(tdata, "genus")` first gives the same embedding. Features
without a genus, or whose genus MGM never saw, are left out with one warning
that counts them; on GlobalPatterns that is 96 of 996 genus-level features.
Each sample then goes through MGM's own preprocessing: relative abundances
over the genera MGM knows, standardised with MGM's per-genus mean and standard
deviation, genera sorted by that value, and the sentence `<bos>`, genera,
`<eos>` cut to 512 tokens. A sample with none of MGM's genera is embedded
from `<bos> <eos>` alone, with a warning naming it.

The embedding is the mean of the last hidden layer over the sample's tokens,
the "element-wise mean pooling" MGM's paper uses for the pretrained model:
256 float32 values per sample, within 2e-6 of MGM 0.5.8's own code. The model
runs on the CPU, one sample at a time; on GlobalPatterns' genus profiles that
is about 60 samples a second on 8 threads, with no memory beyond the model's.
On a CPU running more than four threads, the first call in a session can
differ from later ones by up to about 2e-4: torch 2.13 and 2.14 sometimes
compute their first `tanh` less precisely.

Cite MGM when you publish results that use it: Zhang H, Zhang Y, Kang Z,
Xiong J, Yang R, Ning K (2026) MGM as a large-scale pretrained foundation
model for microbiome analyses in diverse contexts. *Adv Sci* 13:e13333.

### Filtering and leakage

A model reads what it was trained on, so filter first: a later feature change
(`bt.pp.filter_features`, `bt.pp.tax_glom`) drops `obsm`, and the embedding
with it. MGM embeds each sample from that sample alone, so computing its
embedding before a cross-validation split leaks nothing. A plugin that
normalises across samples would leak; check how yours works.

### Adding a model

Models are plugins. A package provides a function that takes the AnnData and
returns a samples x dimensions NumPy array, leaving the AnnData unchanged, and
registers it under a name in the entry-point group `biotapy.embeddings` of its
`pyproject.toml`:

```toml
[project.entry-points."biotapy.embeddings"]
mymodel = "mypackage.embedding:embed"
```

Once the package is installed, `bt.ml.embed(adata, "mymodel")` finds it;
nothing is registered by hand. biotapy checks what the plugin returns - a
finite float array with one row per sample - and raises an error naming the
plugin otherwise. An unknown name raises `KeyError` listing the installed
models, and a name that two installed packages register raises instead of
picking one.
