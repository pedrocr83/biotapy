"""The one download cache: every file biotapy fetches at run time goes through a pooch made here."""

from typing import cast

import pooch


def make_pooch(base_url: str, registry: dict[str, str | None], *, urls: dict[str, str] | None = None) -> pooch.Pooch:
    """A pooch over biotapy's cache: ``BIOTAPY_DATA_DIR`` if set, else pooch's per-user cache directory."""
    # pooch ships no py.typed marker, so mypy --strict infers Any for the untyped `create`; cast it back to Pooch.
    return cast(
        pooch.Pooch,
        pooch.create(
            path=pooch.os_cache("biotapy"), base_url=base_url, registry=registry, urls=urls, env="BIOTAPY_DATA_DIR"
        ),
    )
