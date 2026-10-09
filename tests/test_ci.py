"""CI must run the rules.md gate (prek hooks), not only local pre-commit."""

import ast
import json
import tomllib
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = yaml.safe_load((ROOT / ".github/workflows/test.yaml").read_text(encoding="utf-8"))
DOCS = ROOT / "docs"


def test_ci_runs_every_prek_hook():
    steps = WORKFLOW["jobs"]["lint"]["steps"]
    assert any("prek run --all-files" in step.get("run", "") for step in steps)


def test_rules_gate_blocks_merges():
    assert "lint" in WORKFLOW["jobs"]["check"]["needs"]


def test_no_extras_job_imports_every_submodule():
    steps = WORKFLOW["jobs"]["import-without-extras"]["steps"]
    assert any("walk_packages" in step.get("run", "") for step in steps)


def test_network_job_runs_network_and_golden_tests():
    steps = WORKFLOW["jobs"]["network"]["steps"]
    assert any(step.get("run", "").strip() == 'uv run --group test pytest -m "network or golden"' for step in steps)


def test_network_job_blocks_merges():
    assert "network" in WORKFLOW["jobs"]["check"]["needs"]


def test_r_bridge_job_runs_the_r_marker():
    steps = WORKFLOW["jobs"]["r-bridge"]["steps"]
    run = [step for step in steps if step.get("run", "").strip() == "uv run --group test --extra r pytest -m r"]
    assert run and "BIOTAPY_DATA_DIR" in run[0]["env"] and run[0]["env"]["RPY2_CFFI_MODE"] == "API"


def test_r_bridge_job_installs_the_image_r_and_packages():
    steps = {step.get("uses", "").split("@")[0]: step for step in WORKFLOW["jobs"]["r-bridge"]["steps"]}
    assert steps["r-lib/actions/setup-r"]["with"]["r-version"] == "4.5.3"
    dockerfile = (ROOT / "tests" / "r" / "Dockerfile").read_text(encoding="utf-8")
    assert "FROM rocker/r-ver:4.5.3" in dockerfile and 'version = "3.22"' in dockerfile
    packages = steps["r-lib/actions/setup-r-dependencies"]["with"]["packages"]
    assert {name.strip() for name in packages.split(",")} == {"bioc::ALDEx2", "bioc::maaslin3"}


def test_r_bridge_job_uses_the_image_cran_snapshot_and_has_a_timeout():
    job = WORKFLOW["jobs"]["r-bridge"]
    setup = next(step for step in job["steps"] if step.get("uses", "").startswith("r-lib/actions/setup-r@"))["with"]
    # The snapshot comes from the rocker/r-ver:4.5.3 base image (tests/r/Dockerfile names no URL): `R -e 'getOption("repos")'`.
    assert setup["cran"] == "https://p3m.dev/cran/__linux__/noble/2026-04-23"
    assert setup["use-public-rspm"] is False
    assert job["timeout-minutes"] == 30


def test_r_bridge_job_blocks_merges():
    assert "r-bridge" in WORKFLOW["jobs"]["check"]["needs"]


def test_ml_extras_job_runs_the_torch_marker_on_python_3_13():
    job = WORKFLOW["jobs"]["ml-extras"]
    setup = next(step for step in job["steps"] if step.get("uses", "").startswith("astral-sh/setup-uv@"))
    assert job["runs-on"] == "ubuntu-latest" and setup["with"]["python-version"] == "3.13"
    runs = [step.get("run", "").strip() for step in job["steps"]]
    assert "uv run --group test --extra torch pytest -m torch" in runs


def test_ml_extras_job_blocks_merges():
    assert "ml-extras" in WORKFLOW["jobs"]["check"]["needs"]


def test_the_torch_extra_comes_from_the_cpu_index_and_nothing_else_does():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert pyproject["project"]["optional-dependencies"]["torch"] == ["torch>=2.9"]
    uv = pyproject["tool"]["uv"]
    assert uv["sources"] == {"torch": {"index": "pytorch-cpu"}}
    assert uv["index"] == [{"name": "pytorch-cpu", "url": "https://download.pytorch.org/whl/cpu", "explicit": True}]


def test_coverage_below_90_percent_fails_the_test_job():
    coverage = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["coverage"]
    assert coverage["report"]["fail_under"] == 90
    steps = WORKFLOW["jobs"]["test"]["steps"]
    assert any(":cov-report" in step.get("run", "") for step in steps)


def test_docs_job_builds_the_docs_with_the_pooch_cache():
    steps = WORKFLOW["jobs"]["docs"]["steps"]
    build = [step for step in steps if step.get("run", "").strip() == "uvx hatch run docs:build"]
    assert build and "BIOTAPY_DATA_DIR" in build[0]["env"]


def test_docs_job_blocks_merges():
    assert "docs" in WORKFLOW["jobs"]["check"]["needs"]


def test_lint_job_imports_the_benchmarks():
    steps = WORKFLOW["jobs"]["lint"]["steps"]
    check = [step for step in steps if step.get("run", "").strip() == "uv run --group dev asv check --python=same"]
    assert check and check[0]["working-directory"] == "benchmarks"


def _conf_constants():
    # Parsed, not imported: conf.py imports Sphinx extensions the test group does not install.
    tree = ast.parse((DOCS / "conf.py").read_text(encoding="utf-8"))
    return {
        target.id: ast.literal_eval(node.value)
        for node in tree.body
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant)
        for target in node.targets
        if isinstance(target, ast.Name)
    }


def test_docs_build_executes_notebooks_and_fails_on_a_cell_error():
    conf = _conf_constants()
    assert conf["nb_execution_mode"] == "cache"
    assert conf["nb_execution_raise_on_error"] is True


def test_a_notebook_cell_has_time_to_download_a_dataset():
    assert _conf_constants()["nb_execution_timeout"] == 300


def _sources():
    return [
        path
        for path in sorted(DOCS.rglob("*"))
        if path.suffix in {".md", ".ipynb"} and not {"_build", "generated"} & set(path.relative_to(DOCS).parts)
    ]


def _notebook_metadata(path):
    if path.suffix == ".ipynb":
        return json.loads(path.read_text(encoding="utf-8")).get("metadata", {})
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return {}
    return yaml.safe_load(text[4 : text.index("\n---", 3)]) or {}


def test_no_page_overrides_the_notebook_execution_settings():
    assert _sources()
    overridden = [
        str(path.relative_to(DOCS))
        for path in _sources()
        if "mystnb" in _notebook_metadata(path) or "execution_mode" in str(_notebook_metadata(path))
    ]
    assert overridden == []
