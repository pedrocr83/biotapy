from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.spatial import procrustes

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def _with_braycurtis():
    tdata = bt.datasets.global_patterns()
    bt.tl.beta(tdata, inplace=True)
    return tdata


def test_pcoa_matches_ape_pcoa():
    vectors = pd.read_csv(GOLDEN / "pcoa_braycurtis_vectors.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
    values = pd.read_csv(GOLDEN / "pcoa_braycurtis_values.csv.gz")
    coords, axes = bt.tl.pcoa(_with_braycurtis(), n_components=10)
    assert list(coords.index) == list(vectors.index)
    np.testing.assert_allclose(axes["eigenvalue"], values["eigenvalue"], rtol=1e-6)
    np.testing.assert_allclose(axes["proportion_explained"], values["relative_eig"], rtol=1e-6)
    ours, theirs = coords.to_numpy(), vectors.to_numpy()
    # Eigenvector signs are arbitrary: align each axis before comparing.
    signs = np.sign((ours * theirs).sum(axis=0))
    np.testing.assert_allclose(ours * signs, theirs, rtol=1e-6)


def test_nmds_matches_vegan_metamds():
    points = pd.read_csv(GOLDEN / "nmds_braycurtis_points.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
    r_stress = pd.read_csv(GOLDEN / "nmds_braycurtis_stress.csv.gz")["stress"].iloc[0]
    coords, stress = bt.tl.nmds(_with_braycurtis(), seed=0)
    assert list(coords.index) == list(points.index)
    assert abs(stress - r_stress) < 0.02
    _, _, disparity = procrustes(points.to_numpy(), coords.to_numpy())
    assert np.sqrt(1 - disparity) > 0.95
