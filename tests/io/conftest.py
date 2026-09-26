import json
import uuid
import zipfile

import biom
import numpy as np
import pytest
import scipy.sparse as sp
from biom.util import biom_open

NEWICK = "((OTU_1:0.1,OTU_2:0.2):0.05,(OTU_3:0.3,OTU_4:0.4):0.1);"
TAXONOMY = [
    ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__", "f__", "g__", "s__"],
    ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__", "g__", "s__"],
    ["k__Bacteria", "p__Bacteroidetes", "c__", "o__", "f__", "g__", "s__"],
    ["k__Archaea", "p__Euryarchaeota", "c__", "o__", "f__", "g__", "s__"],
]


@pytest.fixture
def biom_table() -> biom.Table:
    """4 features x 3 samples, laid out as BIOM stores them."""
    counts = sp.csr_matrix(np.array([[5, 0, 3], [1, 2, 0], [0, 4, 6], [7, 0, 0]]))
    return biom.Table(
        counts,
        ["OTU_1", "OTU_2", "OTU_3", "OTU_4"],
        ["S1", "S2", "S3"],
        observation_metadata=[{"taxonomy": t} for t in TAXONOMY],
        sample_metadata=[{"group": g} for g in ["A", "A", "B"]],
    )


@pytest.fixture
def biom_hdf5(tmp_path, biom_table):
    path = tmp_path / "table.biom"
    with biom_open(str(path), "w") as handle:
        biom_table.to_hdf5(handle, "biotapy tests")
    return path


@pytest.fixture
def biom_json(tmp_path, biom_table):
    path = tmp_path / "table.json.biom"
    path.write_text(biom_table.to_json("biotapy tests"))
    return path


@pytest.fixture
def newick(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text(NEWICK)
    return path


@pytest.fixture
def biom_json_ids(tmp_path):
    """Write a minimal BIOM 1.0 JSON table with ids exactly as given (ints allowed)."""

    def write(sample_ids, observation_ids):
        path = tmp_path / "ids.biom"
        document = {
            "id": None,
            "format": "Biological Observation Matrix 1.0.0",
            "format_url": "http://biom-format.org",
            "type": "OTU table",
            "generated_by": "biotapy tests",
            "date": "2026-09-26T00:00:00",
            "rows": [{"id": i, "metadata": None} for i in observation_ids],
            "columns": [{"id": i, "metadata": None} for i in sample_ids],
            "matrix_type": "sparse",
            "matrix_element_type": "int",
            "shape": [len(observation_ids), len(sample_ids)],
            "data": [[0, 0, 5], [1, 1, 3]],
        }
        path.write_text(json.dumps(document))
        return path

    return write


@pytest.fixture
def make_qza(tmp_path):
    """Build a minimal .qza: <uuid>/metadata.yaml, <uuid>/VERSION, <uuid>/data/<payload>."""

    def make(name, payload, content, semantic_type):  # noqa: PLR0917 -- one slot per .qza concept; test call sites are positional
        uid = str(uuid.uuid5(uuid.NAMESPACE_URL, name))  # deterministic (R11.4)
        path = tmp_path / f"{name}.qza"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(f"{uid}/metadata.yaml", f"uuid: {uid}\ntype: {semantic_type}\nformat: null\n")
            archive.writestr(f"{uid}/VERSION", "QIIME 2\narchive: 5\nframework: 2024.10\n")
            archive.writestr(f"{uid}/data/{payload}", content)
        return path

    return make
