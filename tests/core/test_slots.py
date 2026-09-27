import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from biotapy._core import add_provenance, feature_subset, infer_x_kind, require_counts, x_kind


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


def test_require_counts_rejects_relative():
    adata = _adata()
    adata.uns["biotapy"] = {"x_kind": "relative"}
    with pytest.raises(ValueError, match="pp.rarefy needs raw counts"):
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
