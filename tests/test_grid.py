#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of Shareloc
# (see https://github.com/CNES/shareloc).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
Module to test functions that use direct grid
"""


import os
import pytest
import numpy as np

from utils import TEST_DIR
from shareloc.grid import Grid


@pytest.mark.unit_tests
def test_grid_geotiff():
    """
    test grid readers
    """
    geotiff_grid_path = os.path.join(TEST_DIR, "ellipsoide", "loc_direct_grid_PHR_2013072139303958CP.tif")
    gri_geotiff = Grid(geotiff_grid_path)
    res_geotiff = gri_geotiff.direct_loc_h(0.5, 0.5, 1000.0)
    np.testing.assert_allclose(res_geotiff, [2.183908972985368, 48.94317692547565, 1000.0], rtol=0, atol=1e-9)
    res_geotiff = gri_geotiff.direct_loc_h(50, 100, 200.0)
    np.testing.assert_allclose(res_geotiff, [2.1828713504608683, 48.942429997483146, 200.0], rtol=0, atol=1e-9)
