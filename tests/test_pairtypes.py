import pytest
import soluanalysis as solu
from pathlib import Path


def test_pairs():
    """
    Test that pair types work.
    """
    pair13 = solu.james.Pair(1, 3)
    pair31 = solu.james.Pair(3, 1)
    assert pair13 == pair31
