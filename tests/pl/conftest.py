import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure


@pytest.fixture
def ax() -> Axes:
    """Axes on a figure pyplot does not track, so tests leave no open figures behind."""
    return Figure().add_subplot()
