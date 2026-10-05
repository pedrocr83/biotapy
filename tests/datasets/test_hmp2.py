import gzip

import mudata
import numpy as np
import pandas as pd
import pytest

import biotapy as bt
from biotapy.datasets import _hmp2, _remote

# Synthetic files in the layout of the three HMP2 products. No HMP2 data is copied: the IBDMDB states no licence.
_METADATA = pd.DataFrame(
    {
        "Project": "HMP2",
        "External ID": ["S1A_P", "S1B_P", "S1T_P", "S2B", "S2A", "S3A_P"],
        "Participant ID": ["C3001", "C3001", "C3001", "C3002", "C3002", "H4001"],
        "data_type": [
            "metagenomics",
            "metagenomics",
            "metatranscriptomics",
            "metagenomics",
            "metagenomics",
            "metagenomics",
        ],
        "week_num": [2.0, 1.0, 0.0, 4.0, 4.0, 0.0],
        "visit_num": [2, 1, 0, 4, 4, 1],
        "diagnosis": ["UC", "UC", "UC", "CD", "CD", "nonIBD"],
        "site_name": ["Cedars-Sinai", "Cedars-Sinai", "Cedars-Sinai", "Cedars-Sinai", "Cedars-Sinai", "MGH"],
        "sex": ["Female", "Female", "Female", "Male", "Male", "Female"],
        "consent_age": [30.0, 30.0, 30.0, 52.0, 52.0, np.nan],
        "Antibiotics": ["No", "No", "No", "Yes", "Yes", "No"],
        "hbi": [np.nan, 3.0, np.nan, 5.0, 5.0, np.nan],
    }
)
_METAGENOMES = ["S1A_P", "S1B_P", "S2B", "S2A", "S3A_P"]
# CPM rows; PWY-2 is found only in S2B, which hmp2() does not keep.
_PATHWAYS = {
    "UNMAPPED": [300000, 250000, 200000, 400000, 350000],
    "UNINTEGRATED": [600000, 650000, 700000, 500000, 550000],
    "PWY-1": [100000, 100000, 50000, 100000, 100000],
    "PWY-1|g__Anaerostipes.s__Anaerostipes_hadrus": [60000, 70000, 0, 40000, 100000],
    "PWY-1|unclassified": [40000, 30000, 0, 60000, 0],
    "PWY-2": [0, 0, 50000, 0, 0],
}
# MetaPhlAn 3 percentages; the leaves (UNKNOWN and two species) of each sample sum to 100.
_LINEAGE = "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Anaerostipes"
_SPECIES = {
    "s__Anaerostipes_hadrus": [60.0, 50.0, 0.0, 30.0, 90.0],
    "s__Anaerostipes_caccae": [30.0, 45.0, 100.0, 70.0, 0.0],
}
_UNKNOWN = [10.0, 5.0, 0.0, 0.0, 10.0]


def _taxa_rows():
    clade = np.add(*_SPECIES.values())
    parts = _LINEAGE.split("|")
    rows = {"UNKNOWN": _UNKNOWN} | {"|".join(parts[: depth + 1]): clade for depth in range(len(parts))}
    return rows | {f"{_LINEAGE}|{name}": values for name, values in _SPECIES.items()}


def _write_table(path, rows, suffix):
    lines = ["Feature\\Sample\t" + "\t".join(f"{sample}{suffix}" for sample in _METAGENOMES)]
    lines += [f"{row}\t" + "\t".join(str(value) for value in values) for row, values in rows.items()]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


@pytest.fixture
def fetched(tmp_path, monkeypatch):
    # As in test_remote.py and test_enzyme.py: the only offline route to the loader is its private _fetch (R11.4).
    _METADATA.to_csv(tmp_path / "hmp2_metadata_2018-08-20.csv", index=False)
    _write_table(tmp_path / "pathabundances_3.tsv.gz", _PATHWAYS, "_pathabundance_cpm")
    _write_table(tmp_path / "taxonomic_profiles_3.tsv.gz", _taxa_rows(), "_profile")
    names = []
    monkeypatch.setattr(_hmp2, "_fetch", lambda name: names.append(name) or str(tmp_path / name))
    return names


def test_keeps_each_participants_first_metagenome(fetched):
    mdata = bt.datasets.hmp2()
    # S1T_P is a metatranscriptome; S2A and S2B share week 4, so the External ID decides.
    assert mdata.obs_names.tolist() == ["S1B_P", "S2A", "S3A_P"]


def test_the_earlier_visit_wins_a_same_week_tie(fetched, tmp_path):
    # S2A has the lower External ID but the later visit; visit_num is used for ordering only.
    _METADATA.assign(visit_num=[2, 1, 0, 4, 5, 1]).to_csv(tmp_path / "hmp2_metadata_2018-08-20.csv", index=False)
    mdata = bt.datasets.hmp2()
    assert mdata.obs_names.tolist() == ["S1B_P", "S2B", "S3A_P"]
    assert "visit_num" not in mdata.obs.columns


def test_a_missing_visit_number_sorts_last(fetched, tmp_path):
    _METADATA.assign(visit_num=[2, 1, 0, 4, np.nan, 1]).to_csv(tmp_path / "hmp2_metadata_2018-08-20.csv", index=False)
    assert bt.datasets.hmp2().obs_names.tolist() == ["S1B_P", "S2B", "S3A_P"]


def test_three_modalities_over_the_same_samples(fetched):
    mdata = bt.datasets.hmp2()
    assert list(mdata.mod) == ["function", "function_by_taxon", "taxa"]
    assert {key: mod.shape for key, mod in mdata.mod.items()} == {
        "function": (3, 4),
        "function_by_taxon": (3, 2),
        "taxa": (3, 3),
    }
    for mod in mdata.mod.values():
        assert mod.obs_names.tolist() == mdata.obs_names.tolist()


def test_metadata_is_in_the_global_obs_and_every_modality(fetched):
    mdata = bt.datasets.hmp2()
    assert mdata.obs.columns.tolist() == _hmp2.COLUMNS
    assert mdata.obs["diagnosis"].cat.categories.tolist() == ["nonIBD", "UC", "CD"]
    assert mdata.obs["diagnosis"].tolist() == ["UC", "CD", "nonIBD"]
    assert pd.isna(mdata.obs.loc["S3A_P", "consent_age"])
    for mod in mdata.mod.values():
        pd.testing.assert_frame_equal(mod.obs, mdata.obs)


def test_values_and_units_come_from_the_readers(fetched):
    mdata = bt.datasets.hmp2()
    assert mdata["function"].uns["biotapy"]["x_kind"] == "cpm"
    assert mdata["taxa"].uns["biotapy"]["x_kind"] == "relative"
    np.testing.assert_array_equal(mdata["function"][:, "PWY-1"].X.toarray().ravel(), [100000, 100000, 100000])
    np.testing.assert_allclose(np.asarray(mdata["taxa"].X.sum(axis=1)).ravel(), 1.0)
    species = mdata["function_by_taxon"].var["species"]
    assert species.iloc[0] == "Anaerostipes_hadrus" and pd.isna(species.iloc[1])


def test_features_absent_from_the_kept_samples_stay(fetched):
    pwy2 = bt.datasets.hmp2()["function"][:, "PWY-2"]
    assert pwy2.X.nnz == 0


def test_fetches_the_three_pinned_files(fetched):
    bt.datasets.hmp2()
    assert sorted(fetched) == ["hmp2_metadata_2018-08-20.csv", "pathabundances_3.tsv.gz", "taxonomic_profiles_3.tsv.gz"]
    assert all(_remote._REGISTRY[name].startswith("sha256:") for name in fetched)


def test_round_trips_through_h5mu(fetched, tmp_path):
    mdata = bt.datasets.hmp2()
    mdata.write_h5mu(tmp_path / "hmp2.h5mu")
    back = mudata.read_h5mu(tmp_path / "hmp2.h5mu")
    assert list(back.mod) == list(mdata.mod) and back.obs_names.tolist() == mdata.obs_names.tolist()
    assert back.obs["diagnosis"].tolist() == mdata.obs["diagnosis"].tolist()
    assert (back["function_by_taxon"].X != mdata["function_by_taxon"].X).nnz == 0


def test_feeds_renorm_and_contributions(fetched):
    relab = bt.fn.renorm(bt.datasets.hmp2(), "relab")
    assert list(relab.mod) == ["function", "function_by_taxon", "taxa"]
    shares = bt.fn.contributions(relab["function_by_taxon"], "PWY-1")
    np.testing.assert_allclose(shares.loc["S2A"].to_numpy(), [0.04, 0.06])


@pytest.mark.network
def test_hmp2_downloads_and_loads():
    mdata = bt.datasets.hmp2()
    assert {key: mod.shape for key, mod in mdata.mod.items()} == {
        "function": (130, 478),
        "function_by_taxon": (130, 21635),
        "taxa": (130, 579),
    }
    assert mdata.obs["diagnosis"].value_counts().to_dict() == {"CD": 65, "UC": 38, "nonIBD": 27}
    assert mdata.obs["Participant ID"].is_unique
