"""Write the HUMAnN golden files (contracts/r-golden-parity).

Run from the repository root, never in CI:

    uv run --no-project --with humann==3.9 --with pandas==3.0.6 python tests/humann/export_golden.py

HUMAnN's utility scripts need no database for a custom mapping file. Each run
rewrites tests/golden/humann/*.csv.gz (samples as rows, row ids without names)
and tests/golden/humann/VERSIONS.txt.
"""

import subprocess
import tempfile
from importlib.metadata import version
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data" / "humann"
GOLDEN = ROOT / "tests" / "golden" / "humann"
RUNS = {
    "regroup_sum": ("humann_regroup_table", "genefamilies.tsv", ["-c", str(DATA / "regroup_map.tsv")]),
    "regroup_mean": ("humann_regroup_table", "genefamilies.tsv", ["-c", str(DATA / "regroup_map.tsv"), "-f", "mean"]),
    "renorm_relab": ("humann_renorm_table", "pathabundance.tsv", ["-u", "relab"]),
    "renorm_cpm": ("humann_renorm_table", "pathabundance.tsv", ["-u", "cpm"]),
    "renorm_relab_nospecial": ("humann_renorm_table", "pathabundance.tsv", ["-u", "relab", "-s", "n"]),
    "renorm_cpm_genefamilies": ("humann_renorm_table", "genefamilies.tsv", ["-u", "cpm"]),
}


def to_golden(path: Path) -> pd.DataFrame:
    """A HUMAnN table as samples x row ids, names dropped (``ID: name|taxon`` -> ``ID|taxon``)."""
    table = pd.read_csv(path, sep="\t", index_col=0)
    ids = table.index.str.split("|", n=1)
    table.index = [parts[0].split(": ", 1)[0] + ("|" + parts[1] if len(parts) == 2 else "") for parts in ids]
    table.columns = table.columns.str.replace(r"_Abundance(-RPKs)?$", "", regex=True)
    return table.T.rename_axis("sample_id")


def main() -> None:
    GOLDEN.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        for name, (tool, source, options) in RUNS.items():
            out = Path(tmp) / f"{name}.tsv"
            subprocess.run([tool, "-i", str(DATA / source), "-o", str(out), *options], check=True)
            # mtime=0: two runs write bit-identical files (playbooks/regenerate-golden-files).
            to_golden(out).to_csv(GOLDEN / f"{name}.csv.gz", compression={"method": "gzip", "mtime": 0})
    (GOLDEN / "VERSIONS.txt").write_text(f"humann {version('humann')}\npandas {version('pandas')}\n")


if __name__ == "__main__":
    main()
