# Standard imports
import os

# Third party imports
import numpy as np
import pytest
import rasterio

# Shareloc imports
from shareloc.draft_pybind.rpc import RPC


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_phrdimap(col, row, alt):
    """
    test the sequence of a inverse localization followed by a direct localization using dimap file
    """

    file_dimap = "/home/adevaux/Bureau/sharloc_231/shareloc/tests/data/rpc/PHRDIMAP_P1BP--2017030824934340CP.XML"

    fctrat = RPC(file_dimap)

    (lonlatalt) = fctrat.direct_loc_h(row, col, alt)

    (row_ar, col_ar, __) = fctrat.inverse_loc(lonlatalt[0][0], lonlatalt[0][1], lonlatalt[0][2])
    assert col_ar == pytest.approx(col, abs=2e-2)
    assert row_ar == pytest.approx(row, abs=2e-2)
    assert False
