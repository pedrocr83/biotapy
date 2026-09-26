"""CI must run the rules.md gate (prek hooks), not only local pre-commit."""

from pathlib import Path

import yaml

WORKFLOW = yaml.safe_load(
    (Path(__file__).resolve().parents[1] / ".github/workflows/test.yaml").read_text(encoding="utf-8")
)


def test_ci_runs_every_prek_hook():
    steps = WORKFLOW["jobs"]["lint"]["steps"]
    assert any("prek run --all-files" in step.get("run", "") for step in steps)


def test_rules_gate_blocks_merges():
    assert "lint" in WORKFLOW["jobs"]["check"]["needs"]


def test_no_extras_job_imports_every_submodule():
    steps = WORKFLOW["jobs"]["import-without-extras"]["steps"]
    assert any("walk_packages" in step.get("run", "") for step in steps)
