"""scripts/knowledge_stale.sh must fail loudly when it cannot compare against the trunk."""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(sys.platform == "win32", reason="bash script; CI runs it on ubuntu")
def test_touched_mode_fails_on_unknown_ref():
    result = subprocess.run(
        ["bash", "scripts/knowledge_stale.sh", "--touched", "--against", "biotapy-no-such-ref"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 2
    assert "biotapy-no-such-ref" in result.stderr
