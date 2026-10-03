import numpy as np
import pandas as pd
import pytest
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
