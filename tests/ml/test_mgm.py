import re
import subprocess
import sys
from importlib.metadata import entry_points
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors

import biotapy as bt

DATA = Path(__file__).parents[1] / "data" / "mgm"
EMBEDDINGS = Path(__file__).parents[2] / "docs" / "tutorials" / "embeddings.md"
# Phase 4 exit gate 2 on GlobalPatterns' genera. docs/tutorials/embeddings.md quotes these outputs; its build has no torch,
# so an mgm test runs the page's code and a default-run test checks the page quotes them.
GLOBAL_PATTERNS_LEFT_OUT = (
    "mgm leaves out 96 of 996 features: 0 without a genus and 96 whose genus is not one of MGM's "
    "(4-29, 4041AA30, A17, Aquamonas, Arctic95A-2, ...)"
)
GLOBAL_PATTERNS_SHAPE = (26, 256)
# The size of MGM's vocabulary (phylogeny.csv's genera), which the tutorial quotes.
MGM_GENERA = 9665
# torch 2.13-2.14's first tanh in a process can saturate on CPUs running more than 4 threads (measured: the
# embedding's first call off by up to 1.6e-4, later calls by 1.7e-6), so equality is checked to 1e-3.
ATOL = 1e-3
# Only a process's first call is affected by that race, so a later call matches MGM's own output to 1e-5 (measured 1.7e-6);
# a changed last token (an <eos> kept at the cut) moves a row by 7.7e-4, which 1e-3 hides.
WARM_ATOL = 1e-5


def _reference_table():
    """tests/data/mgm/counts.csv as an AnnData: 7 samples x 617 features, var["genus"] (one missing)."""
    counts = pd.read_csv(DATA / "counts.csv", index_col=0)
    return AnnData(
        X=sp.csr_matrix(counts.drop(columns="genus").T.to_numpy()),
        obs=pd.DataFrame(index=counts.columns[1:]),
        var=counts[["genus"]],
    )


def test_mgm_is_registered_as_a_plugin():
    points = [point for point in entry_points(group="biotapy.embeddings") if point.name == "mgm"]
    assert [point.value for point in points] == ["biotapy.ml._mgm:embed"]


def test_without_a_genus_column_raises(make_adata):
    with pytest.raises(KeyError, match=r"mgm reads var\['genus'\], which this table does not have"):
        bt.ml.embed(make_adata(np.array([[1, 2]])), "mgm")


def test_negative_values_raise(make_adata):
    adata = make_adata(np.array([[1.0, -2.0]]))
    adata.var["genus"] = ["Bacteroides", "Prevotella"]
    with pytest.raises(ValueError, match="mgm reads counts or relative abundances; adata.X holds negative"):
        bt.ml.embed(adata, "mgm")


def test_without_the_extra_names_it(monkeypatch):
    # None in sys.modules makes the import fail whether or not the extra is installed.
    monkeypatch.setitem(sys.modules, "torch", None)
    monkeypatch.setitem(sys.modules, "transformers", None)
    with pytest.raises(ImportError, match=r"pip install 'biotapy\[mgm\]'"):
        bt.ml.embed(bt.pp.tax_glom(bt.datasets.toy(), "genus"), "mgm")


@pytest.mark.mgm
def test_matches_mgm_s_own_embedding():
    expected = pd.read_csv(DATA / "embeddings.csv", index_col=0)
    with pytest.warns(UserWarning) as record:
        result = bt.ml.embed(_reference_table(), "mgm")
    assert [str(warning.message) for warning in record] == [
        "mgm leaves out 2 of 617 features: 1 without a genus and 1 whose genus is not one of MGM's (Notagenus)",
        "mgm embeds 2 sample(s) from <bos> <eos> alone, none of their genera being MGM's: s4, s5",
    ]
    assert result.shape == (7, 256) and result.dtype == np.float32
    np.testing.assert_allclose(result, expected.to_numpy(), rtol=0, atol=ATOL)
    with pytest.warns(UserWarning):
        warm = bt.ml.embed(_reference_table(), "mgm")
    np.testing.assert_allclose(warm, expected.to_numpy(), rtol=0, atol=WARM_ATOL)


# A process's first VML call races inside torch's oneMKL (pytorch/pytorch#188792), so each run is a fresh process,
# oversubscribed to 32 threads (about 15% differ without the fix, hence 24 runs); the first sample fills MGM's 512
# tokens, the largest GELU MGM computes.
FIRST_CALL_SCRIPT = """
import sys, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, torch
from anndata import AnnData
import biotapy as bt
torch.set_num_threads(32)
genera = pd.read_csv(sys.argv[1], index_col=0)["genus"].dropna()
counts = np.random.default_rng(0).integers(1, 1000, size=(2, genera.size)).astype(float)
table = AnnData(sp.csr_matrix(counts), obs=pd.DataFrame(index=["a", "b"]),
                var=pd.DataFrame({"genus": list(genera)}, index=[f"f{i}" for i in range(genera.size)]))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    print(int(np.array_equal(bt.ml.embed(table, "mgm"), bt.ml.embed(table, "mgm"))))
"""


@pytest.mark.mgm
def test_a_fresh_process_s_first_embedding_is_its_second():
    for _ in range(24):
        run = subprocess.run(
            [sys.executable, "-c", FIRST_CALL_SCRIPT, str(DATA / "counts.csv")],
            capture_output=True,
            text=True,
            check=True,
        )
        assert run.stdout.strip() == "1"


@pytest.mark.mgm
def test_features_of_one_genus_are_summed_so_tax_glom_changes_nothing():
    tdata = bt.datasets.toy()  # f4 and f5 are both Bacteroides; f8 has no genus
    with pytest.warns(UserWarning, match="1 without a genus"):
        raw = bt.ml.embed(tdata, "mgm")
    np.testing.assert_allclose(raw, bt.ml.embed(bt.pp.tax_glom(tdata, "genus"), "mgm"), rtol=0, atol=ATOL)


@pytest.mark.mgm
def test_the_leave_out_warning_lists_unknown_genera_only_when_there_are_some():
    with pytest.warns(UserWarning) as record:
        bt.ml.embed(bt.datasets.toy(), "mgm")
    assert [str(warning.message) for warning in record] == [
        "mgm leaves out 1 of 8 features: 1 without a genus and 0 whose genus is not one of MGM's"
    ]


@pytest.mark.mgm
def test_relative_abundances_embed_as_their_counts():
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    relative = genera.copy()
    relative.X = sp.csr_matrix(bt.pp.relative(genera).layers["relative"])
    np.testing.assert_allclose(bt.ml.embed(relative, "mgm"), bt.ml.embed(genera, "mgm"), rtol=0, atol=ATOL)


@pytest.mark.mgm
def test_feature_order_does_not_matter():
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    shuffled = genera[:, ::-1].copy()
    np.testing.assert_allclose(bt.ml.embed(shuffled, "mgm"), bt.ml.embed(genera, "mgm"), rtol=0, atol=ATOL)


@pytest.mark.mgm
def test_an_all_zero_feature_changes_nothing():
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    padded = AnnData(
        X=sp.hstack([genera.X, sp.csr_matrix((genera.n_obs, 1), dtype=genera.X.dtype)], format="csr"),
        obs=genera.obs,
        var=pd.DataFrame({"genus": [*genera.var["genus"], "Akkermansia"]}, index=[*genera.var_names, "zero"]),
    )
    np.testing.assert_allclose(bt.ml.embed(padded, "mgm"), bt.ml.embed(genera, "mgm"), rtol=0, atol=ATOL)


@pytest.mark.mgm
def test_single_sample():
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    alone = bt.ml.embed(genera[[3]].copy(), "mgm")
    np.testing.assert_allclose(alone, bt.ml.embed(genera, "mgm")[[3]], rtol=0, atol=ATOL)


@pytest.mark.mgm
def test_keeps_the_input(assert_unchanged):
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    before = genera.copy()
    bt.ml.embed(genera, "mgm")
    assert_unchanged(before, genera)


@pytest.mark.mgm
def test_embeds_global_patterns_end_to_end():
    # Phase 4 exit gate 2: one foundation model plugged in end to end (decisions 10, 17).
    tdata = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
    with pytest.warns(UserWarning, match=re.escape(GLOBAL_PATTERNS_LEFT_OUT)) as record:
        bt.ml.embed(tdata, "mgm", inplace=True)
    assert not [warning for warning in record if "<bos> <eos>" in str(warning.message)]  # no sample is empty
    embedding = tdata.obsm["X_mgm"]
    assert embedding.shape == GLOBAL_PATTERNS_SHAPE and embedding.dtype == np.float32 and np.isfinite(embedding).all()
    with pytest.warns(UserWarning, match="mgm leaves out 96"):
        np.testing.assert_allclose(bt.ml.embed(tdata, "mgm"), embedding, rtol=0, atol=ATOL)
    # Every sample's nearest neighbour in the embedding comes from the same environment.
    neighbour = NearestNeighbors(n_neighbors=2, metric="cosine").fit(embedding).kneighbors(embedding)[1][:, 1]
    types = tdata.obs["SampleType"].to_numpy()
    assert (types[neighbour] == types).all()


@pytest.mark.mgm
def test_the_embedding_tutorial_gives_the_outputs_it_quotes(run_page):
    with pytest.warns(UserWarning, match=re.escape(GLOBAL_PATTERNS_LEFT_OUT)):
        namespace = run_page(EMBEDDINGS)
    assert (namespace["embedding"].shape, namespace["embedding"].dtype) == (GLOBAL_PATTERNS_SHAPE, np.float32)
    assert namespace["same_type"] == GLOBAL_PATTERNS_SHAPE[0]


def test_the_embedding_tutorial_quotes_the_end_to_end_outputs():
    page = EMBEDDINGS.read_text(encoding="utf-8")
    assert f"```text\nUserWarning: {GLOBAL_PATTERNS_LEFT_OUT}\n```" in page
    assert f"```text\n({GLOBAL_PATTERNS_SHAPE}, dtype('float32'))\n```" in page
    assert f"```text\n{GLOBAL_PATTERNS_SHAPE[0]}\n```" in page
    assert (
        f"For all {GLOBAL_PATTERNS_SHAPE[0]} samples, the nearest neighbour is a sample of the same type"
        in " ".join(page.split())
    )


def _cached_file(name):
    """The file ``name`` that MGM's wheel was extracted to, in biotapy's data cache."""
    from biotapy.ml import _mgm

    cache = _mgm.make_pooch(_mgm._BASE_URL, {_mgm._WHEEL: _mgm._SHA256})
    paths = cache.fetch(_mgm._WHEEL, processor=_mgm.pooch.Unzip(members=_mgm._MEMBERS))
    return next(Path(path) for path in paths if Path(path).name == name)


@pytest.mark.mgm
def test_the_embedding_tutorial_s_vocabulary_and_example_genera_hold():
    vocabulary = pd.read_csv(_cached_file("phylogeny.csv"), index_col=0).index
    genera = set(bt.pp.tax_glom(bt.datasets.global_patterns(), "genus").var["genus"].dropna())
    assert len(vocabulary) == MGM_GENERA
    # The clone name and the one-word Candidatus genus the page names: biotapy reads "BD2-13" as MGM's "BD2".
    assert {"BD2-13", "CandidatusPelagibacter"} <= genera
    assert "g__BD2" not in vocabulary and "g__CandidatusPelagibacter" not in vocabulary
    assert "g__Candidatus_Pelagibacter" in vocabulary


def test_the_embedding_tutorial_quotes_its_vocabulary_and_labels_its_timings():
    page = " ".join(EMBEDDINGS.read_text(encoding="utf-8").split())
    samples, dimensions = GLOBAL_PATTERNS_SHAPE
    assert f"MGM's vocabulary holds {MGM_GENERA:,} genera" in page
    assert f"each sample comes back as a vector of {dimensions} numbers" in page
    assert f"{samples} samples" in page
    left_out, features = re.match(r"mgm leaves out (\d+) of (\d+) features", GLOBAL_PATTERNS_LEFT_OUT).groups()
    assert f"{left_out} of GlobalPatterns' {features} genus-level features name a genus outside it" in page
    assert "one unpinned laptop run, not checked by CI" in page


@pytest.mark.mgm
@pytest.mark.parametrize("name", ["phylogeny.csv", "config.json"])
def test_a_damaged_extracted_file_is_extracted_again(name):
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    expected = bt.ml.embed(genera, "mgm")
    path = _cached_file(name)
    original = path.read_bytes()
    path.write_bytes(original[: len(original) // 2])
    result = bt.ml.embed(genera, "mgm")
    assert path.read_bytes() == original
    np.testing.assert_allclose(result, expected, rtol=0, atol=ATOL)


@pytest.mark.mgm
def test_a_file_that_stays_damaged_raises(monkeypatch):
    from biotapy.ml import _mgm

    monkeypatch.setitem(_mgm._SHA256_OF, "phylogeny.csv", "0" * 64)
    genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    with pytest.raises(
        ValueError, match=r"after extracting it again; delete .*microformer_mgm-0\.5\.8-py3-none-any\.whl\.unzip\n"
    ):
        bt.ml.embed(genera, "mgm")
