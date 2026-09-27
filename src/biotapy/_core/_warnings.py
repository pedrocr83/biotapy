"""User-facing warnings attributed to the caller's code, not to biotapy."""

import warnings
from pathlib import Path

# Warnings point at the first frame outside biotapy, however deep the call.
_PACKAGE_DIR = str(Path(__file__).resolve().parents[1])


def warn_user(message: str) -> None:
    """Warn with ``UserWarning`` at the first stack frame outside biotapy."""
    warnings.warn(message, UserWarning, skip_file_prefixes=(_PACKAGE_DIR,))
