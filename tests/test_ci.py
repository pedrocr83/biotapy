"""CI must run the rules.md gate (prek hooks), not only local pre-commit."""

import tomllib
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = yaml.safe_load((ROOT / ".github/workflows/test.yaml").read_text(encoding="utf-8"))


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


def test_coverage_below_90_percent_fails_the_test_job():
    coverage = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["coverage"]
    assert coverage["report"]["fail_under"] == 90
    steps = WORKFLOW["jobs"]["test"]["steps"]
    assert any(":cov-report" in step.get("run", "") for step in steps)
