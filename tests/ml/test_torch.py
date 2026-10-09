import pickle
import sys
import tracemalloc

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt


@pytest.fixture
def torch():
    # Imported here, not at module level: without the extra, collecting this file must still work so the
    # tests below that need no torch run in every job.
    import torch

    return torch


def _stacked(dataset, torch):
    return torch.stack([dataset[i] for i in range(len(dataset))]).numpy()


@pytest.mark.torch
def test_batches_hold_float32_rows_and_int64_labels(torch):
    loader = torch.utils.data.DataLoader(bt.ml.to_torch(bt.datasets.toy(), label_key="group"), batch_size=4)
    features, labels = next(iter(loader))
    assert features.shape == (4, 8) and features.dtype == torch.float32
    assert labels.shape == (4,) and labels.dtype == torch.int64
    assert [len(batch[0]) for batch in loader] == [4, 2]


@pytest.mark.torch
def test_is_a_torch_dataset_with_one_item_per_sample(torch):
    dataset = bt.ml.to_torch(bt.datasets.toy())
    assert isinstance(dataset, torch.utils.data.Dataset) and len(dataset) == 6


@pytest.mark.torch
def test_items_are_the_rows_of_x(torch):
    tdata = bt.datasets.toy()
    np.testing.assert_array_equal(_stacked(bt.ml.to_torch(tdata), torch), tdata.X.toarray().astype(np.float32))


@pytest.mark.torch
def test_labels_follow_the_category_order(torch):
    tdata = bt.datasets.toy()
    tdata.obs["group"] = tdata.obs["group"].cat.reorder_categories(["B", "A"])
    dataset = bt.ml.to_torch(tdata, label_key="group")
    assert [int(dataset[i][1]) for i in range(6)] == [1, 1, 1, 0, 0, 0]


@pytest.mark.torch
def test_a_subset_of_one_dataset_keeps_the_codes_that_per_split_datasets_shift(torch):
    tdata = bt.datasets.toy()
    test = [3, 5]  # both group B: the split lacks class A
    whole = bt.ml.to_torch(tdata, label_key="group")
    subset = torch.utils.data.Subset(whole, test)
    assert [int(subset[i][1]) for i in range(2)] == [1, 1]
    # anndata drops the unused category on subsetting, so a dataset built per split recodes B as 0.
    per_split = bt.ml.to_torch(tdata[test], label_key="group")
    assert [int(per_split[i][1]) for i in range(2)] == [0, 0]


@pytest.mark.torch
@pytest.mark.parametrize(
    ("values", "codes"), [(["b", "a", "c", "a", "b", "c"], [1, 0, 2, 0, 1, 2]), ([True, False] * 3, [1, 0] * 3)]
)
def test_string_and_bool_labels_become_sorted_codes(torch, values, codes):
    tdata = bt.datasets.toy()
    tdata.obs["label"] = values
    dataset = bt.ml.to_torch(tdata, label_key="label")
    assert [int(dataset[i][1]) for i in range(6)] == codes and dataset[0][1].dtype == torch.int64


@pytest.mark.torch
def test_numeric_labels_become_float32(torch):
    tdata = bt.datasets.toy()
    tdata.obs["age"] = [30, 41, 25, 60, 52, 47]
    dataset = bt.ml.to_torch(tdata, label_key="age")
    assert dataset[1][1].dtype == torch.float32 and float(dataset[1][1]) == 41.0


@pytest.mark.torch
@pytest.mark.parametrize(("make", "layer"), [(bt.pp.clr, "clr"), (bt.pp.relative, "relative")])
def test_layer_is_read_instead_of_x(torch, make, layer):
    adata = make(bt.datasets.toy())
    expected = sp.csr_matrix(adata.layers[layer]).toarray().astype(np.float32)
    np.testing.assert_array_equal(_stacked(bt.ml.to_torch(adata, layer=layer), torch), expected)


@pytest.mark.torch
def test_a_view_of_some_samples_gives_those_samples(torch):
    tdata = bt.datasets.toy()
    dataset = bt.ml.to_torch(tdata[[0, 2, 4]], label_key="group")
    rows = torch.stack([dataset[i][0] for i in range(3)]).numpy()
    np.testing.assert_array_equal(rows, tdata.X[[0, 2, 4]].toarray())
    assert [int(dataset[i][1]) for i in range(3)] == [0, 0, 1]


@pytest.mark.torch
def test_keeps_the_input(torch, assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    dataset = bt.ml.to_torch(tdata, label_key="group")
    for features, _ in torch.utils.data.DataLoader(dataset, batch_size=4):
        features += 1
    dataset[0][0][:] = 99
    assert_unchanged(before, tdata)


@pytest.mark.torch
def test_an_item_does_not_share_memory_with_a_dense_float32_layer(torch):
    adata = bt.pp.clr(bt.datasets.toy())
    adata.layers["clr"] = adata.layers["clr"].astype(np.float32)
    before = adata.layers["clr"].copy()
    bt.ml.to_torch(adata, layer="clr")[0][:] = 99
    np.testing.assert_array_equal(adata.layers["clr"], before)


@pytest.mark.torch
def test_the_dataset_survives_pickling_and_spawned_workers(torch):
    # spawn and forkserver (macOS, Windows, Linux on Python 3.14) pickle the dataset to send it to each worker.
    dataset = bt.ml.to_torch(bt.datasets.toy(), label_key="group")
    clone = pickle.loads(pickle.dumps(dataset))
    assert len(clone) == 6 and torch.equal(clone[4][0], dataset[4][0]) and int(clone[4][1]) == int(dataset[4][1])
    workers = torch.utils.data.DataLoader(dataset, batch_size=2, num_workers=2, multiprocessing_context="spawn")
    alone = torch.utils.data.DataLoader(dataset, batch_size=2)
    for (rows, labels), (expected_rows, expected_labels) in zip(workers, alone, strict=True):
        assert torch.equal(rows, expected_rows) and torch.equal(labels, expected_labels)


@pytest.mark.torch
def test_editing_a_returned_label_does_not_change_the_next_fetch(torch):
    dataset = bt.ml.to_torch(bt.datasets.toy(), label_key="group")
    dataset[0][1].add_(3)
    assert int(dataset[0][1]) == 0


@pytest.mark.torch
@pytest.mark.parametrize("sparse", [True, False])
def test_a_slice_index_raises_and_integer_indices_work(torch, sparse):
    tdata = bt.datasets.toy()
    adata = AnnData(X=tdata.X if sparse else tdata.X.toarray())
    dataset = bt.ml.to_torch(adata)
    with pytest.raises(TypeError, match="slice"):
        dataset[1:3]
    expected = tdata.X.toarray().astype(np.float32)
    np.testing.assert_array_equal(dataset[np.int64(2)].numpy(), expected[2])
    np.testing.assert_array_equal(dataset[-1].numpy(), expected[-1])
    with pytest.raises(IndexError):
        dataset[6]


@pytest.mark.torch
def test_references_x_so_a_later_change_shows(torch):
    tdata = bt.datasets.toy()
    dataset = bt.ml.to_torch(tdata)
    tdata.X.data[:] = 0
    assert float(dataset[0].sum()) == 0.0


@pytest.mark.torch
def test_a_large_sparse_table_is_never_dense(torch):
    # 2,000 x 50,000 at 0.1% density: 100,000 stored values; dense it would be 800 MB as float64.
    adata = AnnData(X=sp.random(2_000, 50_000, density=0.001, format="csr", random_state=0))
    tracemalloc.start()
    try:
        dataset = bt.ml.to_torch(adata)
        batch = next(iter(torch.utils.data.DataLoader(dataset, batch_size=64)))
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert batch.shape == (64, 50_000)
    assert peak < 50_000_000


@pytest.mark.torch
def test_all_zero_sample_and_feature_are_zeros(torch, make_adata):
    adata = make_adata(np.array([[0, 0, 0], [3, 0, 1]]))
    rows = _stacked(bt.ml.to_torch(adata), torch)
    np.testing.assert_array_equal(rows, [[0, 0, 0], [3, 0, 1]])


@pytest.mark.torch
def test_single_sample(torch, make_adata):
    dataset = bt.ml.to_torch(make_adata(np.array([[2, 0, 5]])))
    batch = next(iter(torch.utils.data.DataLoader(dataset, batch_size=4)))
    assert len(dataset) == 1 and batch.tolist() == [[2.0, 0.0, 5.0]]


@pytest.mark.torch
def test_other_sparse_formats_are_read_as_csr(torch, make_adata):
    adata = make_adata(np.array([[1, 0], [0, 4]]))
    adata.X = sp.csc_matrix(adata.X)
    np.testing.assert_array_equal(_stacked(bt.ml.to_torch(adata), torch), [[1, 0], [0, 4]])


@pytest.mark.torch
@given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 1000)))
def test_items_stack_back_to_the_table(dense):
    # Hypothesis rejects function-scoped fixtures, so this test imports torch itself.
    import torch

    rows = torch.stack(list(bt.ml.to_torch(AnnData(X=sp.csr_matrix(dense))))).numpy()
    np.testing.assert_array_equal(rows, dense.astype(np.float32))


def test_missing_label_column_raises():
    with pytest.raises(KeyError, match="label_key='diet' is not a column of obs"):
        bt.ml.to_torch(bt.datasets.toy(), label_key="diet")


def test_missing_label_raises():
    tdata = bt.datasets.toy()
    tdata.obs.loc["s2", "group"] = np.nan
    with pytest.raises(ValueError, match=r"label_key='group' has 1 missing value\(s\)"):
        bt.ml.to_torch(tdata, label_key="group")


def test_missing_layer_raises():
    with pytest.raises(KeyError, match=r"layer='clr' is not in adata.layers"):
        bt.ml.to_torch(bt.datasets.toy(), layer="clr")


def test_a_table_that_is_not_an_array_raises():
    adata = AnnData(obs=pd.DataFrame(index=["s1", "s2"]))
    with pytest.raises(TypeError, match="adata.X is a NoneType; to_torch reads a NumPy array or a SciPy sparse matrix"):
        bt.ml.to_torch(adata)


def test_options_are_keyword_only():
    with pytest.raises(TypeError):
        bt.ml.to_torch(bt.datasets.toy(), "group")


def test_without_torch_names_the_extra(monkeypatch):
    # None in sys.modules makes `import torch` fail whether or not the extra is installed.
    monkeypatch.setitem(sys.modules, "torch", None)
    with pytest.raises(ImportError, match=r"pip install 'biotapy\[torch\]'"):
        bt.ml.to_torch(bt.datasets.toy())
