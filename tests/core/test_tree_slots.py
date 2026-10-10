import h5py
import numpy as np
import pytest
from anndata import AnnData

import biotapy as bt
from biotapy._core import TreeData, get_tree, read_tree_slots, write_tree_slots


def _plain(tdata: TreeData) -> AnnData:
    return AnnData(tdata.X.copy(), obs=tdata.obs.copy(), var=tdata.var.copy())


def test_slots_round_trip_through_an_h5_group():
    toy = bt.datasets.toy()
    with h5py.File("slots.h5", "w", driver="core", backing_store=False) as handle:
        group = handle.create_group("mod")
        write_tree_slots(group, toy)
        assert "biotapy-treedata-encoding" in group.attrs
        back = read_tree_slots(group, _plain(toy))
    assert type(back) is TreeData
    assert set(get_tree(back).edges) == set(get_tree(toy).edges)
    np.testing.assert_array_equal(back.X.toarray(), toy.X.toarray())


def test_writing_slots_leaves_the_group_without_other_keys():
    toy = bt.datasets.toy()
    with h5py.File("slots2.h5", "w", driver="core", backing_store=False) as handle:
        group = handle.create_group("mod")
        write_tree_slots(group, toy)
        assert {"obs", "var", "X"}.isdisjoint(group.keys())


def test_an_unmarked_group_raises():
    toy = bt.datasets.toy()
    with h5py.File("slots3.h5", "w", driver="core", backing_store=False) as handle:
        with pytest.raises(KeyError, match="biotapy-treedata-encoding"):
            read_tree_slots(handle.create_group("mod"), _plain(toy))
