# Transforms

Transforms add a same-shape view of your data as a new layer; they never touch `X` or drop
anything else, so you can always get back to what you started from.

## Relative abundance

`bt.pp.relative` scales every sample (row) so its features sum to one, and stores the result in
`layers["relative"]`:

```python
import biotapy as bt

tdata = bt.datasets.toy()
out = bt.pp.relative(tdata)
out.layers["relative"][0].sum()  # 1.0
```

`X` still holds the original counts; `out.uns["biotapy"]["provenance"]` records that
`pp.relative` ran.

### All-zero samples

A sample with no reads has nothing to divide by. Where phyloseq's
`transform_sample_counts` divides by zero and returns `NaN`, biotapy leaves an
all-zero sample as all zeros, so downstream steps do not have to special-case `NaN`.
