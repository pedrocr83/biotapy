import pytest

from biotapy._core import warn_user


def test_warn_user_points_at_the_caller():
    with pytest.warns(UserWarning, match="look here") as record:
        warn_user("look here")
    assert record[0].filename == __file__
