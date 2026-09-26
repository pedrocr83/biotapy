"""Self-test for the `assert_unchanged` purity fixture (rules.md R3.3; F3)."""

import pytest

import biotapy as bt


def test_passes_for_an_untouched_copy(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    assert_unchanged(before, tdata)


def test_catches_provenance_appended_on_the_after_object(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    tdata.uns["biotapy"]["provenance"] = [*tdata.uns["biotapy"]["provenance"], "bogus"]
    with pytest.raises(AssertionError):
        assert_unchanged(before, tdata)


def test_catches_x_kind_change(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    tdata.uns["biotapy"]["x_kind"] = "relative"
    with pytest.raises(AssertionError):
        assert_unchanged(before, tdata)


def test_catches_in_place_layer_edit(assert_unchanged):
    tdata = bt.pp.relative(bt.datasets.toy())
    before = tdata.copy()
    tdata.layers["relative"].data[:] = 0
    with pytest.raises(AssertionError):
        assert_unchanged(before, tdata)
