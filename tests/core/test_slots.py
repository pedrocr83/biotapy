import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from biotapy._core import add_provenance, feature_subset, require_counts, x_kind


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
