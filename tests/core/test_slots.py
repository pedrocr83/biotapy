import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from biotapy._core import add_provenance, feature_subset, infer_x_kind, require_categorical, require_counts, x_kind


def _adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=sp.csr_matrix(np.arange(6).reshape(2, 3)),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["f1", "f2", "f3"]),
    )
    adata.layers["relative"] = adata.X.copy()
    adata.obsm["X_pcoa"] = np.zeros((2, 2))
    adata.obsp["braycurtis"] = sp.csr_matrix((2, 2))
    adata.uns["other"] = 1
    return adata


def test_x_kind_defaults_to_counts():
    assert x_kind(_adata()) == "counts"


def test_require_categorical_rejects_numeric_naming_the_argument():
    values = pd.Series([1.0, 2.0, np.nan], name="dose")
    with pytest.raises(
        TypeError,
        match=r"grouping='dose' is a numeric column \(float64\); tl.permanova compares groups, "
        r"so convert it with \.astype\(\"category\"\) for one group per value",
    ):
        require_categorical(values, arg="grouping", purpose="tl.permanova compares groups")


@pytest.mark.parametrize(
    "values",
    [
        pd.Series([True, False]),
        pd.Series(pd.Categorical([1, 2])),
        pd.Series(["a", None]),
        pd.Series(["a", "b"], dtype="str"),
    ],
)
def test_require_categorical_accepts_bool_category_and_string(values):
    require_categorical(values, arg="x", purpose="plots group by category")


def test_require_counts_rejects_relative():
    adata = _adata()
    adata.uns["biotapy"] = {"x_kind": "relative"}
    with pytest.raises(ValueError, match="pp.rarefy needs raw counts"):
        require_counts(adata, func="pp.rarefy")


def test_require_counts_rejects_non_integer_values():
    adata = _adata()
    adata.X = sp.csr_matrix(np.array([[0.0, 2.0, 1.0], [3.0, 4.0, 5.0]]))
    require_counts(adata, func="pp.rarefy")  # whole numbers stored as float pass
    adata.X = sp.csr_matrix(np.array([[0.0, 0.5, 1.0], [3.0, 4.0, 5.0]]))
    with pytest.raises(ValueError, match="pp.rarefy needs raw counts in X, but X holds non-integer or missing"):
        require_counts(adata, func="pp.rarefy")


def test_require_counts_rejects_nan_values():
    adata = _adata()
    adata.X = sp.csr_matrix(np.array([[0.0, np.nan, 1.0], [3.0, 4.0, 5.0]]))
    with pytest.raises(ValueError, match=r"missing \(NaN\)"):
        require_counts(adata, func="pp.rarefy")


def test_add_provenance_appends_json_entries():
    adata = _adata()
    add_provenance(adata, "pp.a", rank="genus")
    add_provenance(adata, "pp.b")
    entries = [json.loads(e) for e in adata.uns["biotapy"]["provenance"]]
    assert [e["step"] for e in entries] == ["pp.a", "pp.b"]
    assert entries[0]["params"] == {"rank": "genus"}


def test_feature_subset_drops_derived_slots():
    out = feature_subset(_adata(), np.array([0, 2]))
    assert list(out.var_names) == ["f1", "f3"]
    assert not out.layers and not out.obsm and not out.obsp
    assert set(out.uns) == {"biotapy"}


def test_feature_subset_leaves_input_alone():
    adata = _adata()
    feature_subset(adata, np.array([0]))
    assert adata.n_vars == 3 and "relative" in adata.layers and "other" in adata.uns


def test_feature_subset_keeps_x():
    np.testing.assert_array_equal(feature_subset(_adata(), np.array([0, 2])).X.toarray(), [[0, 2], [3, 5]])


def test_infer_x_kind_whole_numbers_are_counts():
    assert infer_x_kind(sp.csr_matrix(np.array([[1.0, 0.0], [2.0, 3.0]]))) == "counts"


def test_infer_x_kind_reads_dense_input():
    assert infer_x_kind(np.array([[4, 0], [0, 7]])) == "counts"


def test_infer_x_kind_rows_summing_to_one_are_relative():
    # the zero row is ignored; 0.99986 is an enterotype-style rounded row sum
    X = sp.csr_matrix(np.array([[0.25, 0.75], [0.0, 0.0], [0.49986, 0.5], [0.0005, 1.0]]))
    assert infer_x_kind(X) == "relative"


def test_infer_x_kind_other_values_are_abundance():
    assert infer_x_kind(sp.csr_matrix(np.array([[0.5, 2.25], [1.5, 0.0]]))) == "abundance"


def test_infer_x_kind_row_sum_outside_tolerance_is_abundance():
    assert infer_x_kind(sp.csr_matrix(np.array([[0.25, 0.75], [0.4, 0.598]]))) == "abundance"


def test_feature_subset_drops_ordination_metadata():
    adata = _adata()
    adata.uns["biotapy"] = {
        "x_kind": "counts",
        "provenance": ["{}"],
        "pcoa": {"eigenvalues": np.ones(2)},
        "nmds": {"stress": 0.1},
    }
    assert feature_subset(adata, np.array([0])).uns["biotapy"] == {"x_kind": "counts", "provenance": ["{}"]}
