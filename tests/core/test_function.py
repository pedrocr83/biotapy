import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from mudata import MuData

from biotapy._core import BY_TAXON_KEY, FUNCTION_KEY, function_var, make_function_mudata

IDS = pd.Index(
    [
        "UNMAPPED",
        "PWY-1: first pathway",
        "PWY-1: first pathway|g__Bacteroides.s__Bacteroides_ovatus",
        "PWY-1: first pathway|unclassified",
        "K1|g__Blautia.s__Blautia_obeum.t__SGB4810",
        "UNINTEGRATED|bug1",
    ]
)


def test_function_var_splits_id_name_and_taxon():
    var = function_var(IDS)
    assert var.index.tolist() == [
        "UNMAPPED",
        "PWY-1",
        "PWY-1|g__Bacteroides.s__Bacteroides_ovatus",
        "PWY-1|unclassified",
        "K1|g__Blautia.s__Blautia_obeum.t__SGB4810",
        "UNINTEGRATED|bug1",
    ]
    assert var["function"].tolist() == ["UNMAPPED", "PWY-1", "PWY-1", "PWY-1", "K1", "UNINTEGRATED"]
    assert var["name"].tolist()[1:4] == ["first pathway"] * 3 and var["name"].isna().tolist()[4:] == [True, True]


def test_function_var_parses_genus_and_species_from_the_stratum():
    var = function_var(IDS)
    assert var["genus"].tolist()[2] == "Bacteroides" and var["species"].tolist()[2] == "Bacteroides_ovatus"
    # HUMAnN 4 strata end in .t__SGB<id>; it is not a species name.
    assert var["species"].tolist()[4] == "Blautia_obeum"
    # "unclassified" and a bare label are a taxon with no rank.
    assert var[["genus", "species"]].iloc[[0, 1, 3, 5]].isna().all().all()


def test_function_var_flags_specials():
    assert function_var(IDS)["special"].tolist() == [True, False, False, False, False, True]


def test_function_var_text_columns_are_the_str_dtype():
    var = function_var(IDS)
    assert all(isinstance(var[column].dtype, pd.StringDtype) for column in ["function", "name", "taxon", "genus"])


def test_function_var_rejects_two_bars():
    with pytest.raises(ValueError, match=r"one '\|'"):
        function_var(pd.Index(["K1|g__A|extra"]))


def test_make_function_mudata_splits_community_and_strata():
    X = np.arange(12, dtype=np.float64).reshape(2, 6)
    mdata = make_function_mudata(X, obs=pd.DataFrame(index=["s1", "s2"]), row_ids=IDS, x_kind="cpm", source="test")
    assert isinstance(mdata, MuData)
    assert mdata[FUNCTION_KEY].var_names.tolist() == ["UNMAPPED", "PWY-1"]
    assert list(mdata[FUNCTION_KEY].var.columns) == ["name", "special"]
    assert list(mdata[BY_TAXON_KEY].var.columns) == ["function", "name", "taxon", "genus", "species", "special"]
    np.testing.assert_array_equal(mdata[BY_TAXON_KEY].X.toarray(), X[:, 2:])
    assert mdata[FUNCTION_KEY].uns["biotapy"]["x_kind"] == "cpm"


def test_make_function_mudata_rejects_repeated_samples():
    with pytest.raises(ValueError, match="duplicate obs ids"):
        make_function_mudata(
            np.ones((2, 1)), obs=pd.DataFrame(index=["s1", "s1"]), row_ids=pd.Index(["K1"]), x_kind="rpk", source="t"
        )


def test_make_function_mudata_round_trips_through_h5mu(tmp_path):
    import mudata

    mdata = make_function_mudata(
        np.ones((2, 6)), obs=pd.DataFrame(index=["s1", "s2"]), row_ids=IDS, x_kind="rpk", source="test"
    )
    mdata.write_h5mu(tmp_path / "f.h5mu")
    back = mudata.read_h5mu(tmp_path / "f.h5mu")
    assert back[BY_TAXON_KEY].var["taxon"].tolist() == mdata[BY_TAXON_KEY].var["taxon"].tolist()
    assert back[FUNCTION_KEY].var["name"].isna().tolist() == [True, False]


def _split(row_ids):
    return make_function_mudata(
        np.ones((2, len(row_ids))),
        obs=pd.DataFrame(index=["s1", "s2"]),
        row_ids=pd.Index(row_ids),
        x_kind="cpm",
        source="test",
    )


def test_make_function_mudata_without_strata_gives_an_empty_by_taxon_modality():
    mdata = _split(["K1", "K2"])
    assert mdata[BY_TAXON_KEY].shape == (2, 0)
    assert mdata[FUNCTION_KEY].shape == (2, 2)


def test_make_function_mudata_with_only_strata_gives_an_empty_function_modality():
    mdata = _split(["K1|g__A.s__B", "K2|unclassified"])
    assert mdata[FUNCTION_KEY].shape == (2, 0)
    assert mdata[BY_TAXON_KEY].shape == (2, 2)


def test_make_function_mudata_global_obs_has_no_columns():
    assert _split(IDS.tolist()).obs.shape[1] == 0


def test_make_function_mudata_stores_x_as_csr_in_both_modalities():
    mdata = _split(IDS.tolist())
    assert all(isinstance(mdata[key].X, sp.csr_matrix) for key in (FUNCTION_KEY, BY_TAXON_KEY))


def test_make_function_mudata_records_x_kind_and_one_provenance_entry_per_modality():
    mdata = _split(IDS.tolist())
    for key in (FUNCTION_KEY, BY_TAXON_KEY):
        meta = mdata[key].uns["biotapy"]
        assert meta["x_kind"] == "cpm"
        assert len(meta["provenance"]) == 1 and '"step": "test"' in meta["provenance"][0]


def test_make_function_mudata_var_names_are_unique_across_modalities():
    names = _split(IDS.tolist()).var_names
    assert names.is_unique and len(names) == len(IDS)
