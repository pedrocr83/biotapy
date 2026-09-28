import subprocess
import sys


def test_import_biotapy_does_not_load_matplotlib():
    # pl imports matplotlib inside its functions, so `import biotapy` stays as fast as before (Task 1.18).
    code = "import sys, biotapy; print('matplotlib' in sys.modules)"
    assert (
        subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip()
        == "False"
    )
