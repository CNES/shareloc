#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2022 Centre National d'Etudes Spatiales (CNES).
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

# Third party imports
import numpy as np
import pytest
import scipy

# Shareloc imports
from shareloc.geomodels.grid import Grid
from shareloc.image import Image

# Shareloc test imports
from tests.helpers import data_path


@pytest.fixture(name="get_geotiff_grid")
def fixture_get_geotiff_grid():
    """
    get grid and associated path
    """
    geotiff_grid_path = data_path("ellipsoide", "loc_direct_grid_PHR_2013072139303958CP.tif")
    gri_geotiff = Grid(geotiff_grid_path)
    grid_image = Image(geotiff_grid_path, read_data=True)

    return gri_geotiff, grid_image


@pytest.mark.unit_tests
def test_grid_geotiff(get_geotiff_grid):
    """
    test grid readers
    """
    gri_geotiff, _ = get_geotiff_grid
    res_geotiff = gri_geotiff.direct_loc_h(0.5, 0.5, 1000.0)
    np.testing.assert_allclose(res_geotiff, [[2.183908972985368, 48.94317692547565, 1000.0]], rtol=0, atol=1e-9)
    res_geotiff = gri_geotiff.direct_loc_h(50, 100, 200.0)
    np.testing.assert_allclose(res_geotiff, [[2.1828713504608683, 48.942429997483146, 200.0]], rtol=0, atol=1e-9)


@pytest.mark.unit_tests
def test_grid_extrapolation(get_geotiff_grid):
    """
    Compare interp with scipy
    """
    # Load grid with Shareloc
    gri_geotiff, grid_image = get_geotiff_grid

    # Load grid with rasterio at 1000m altitude
    data_lon_1000 = grid_image.data[5, :, :]
    # Grid parameters
    grid_step = np.array([grid_image.pixel_size_row, grid_image.pixel_size_col])
    grid_origin = np.array(
        [
            grid_image.origin_row,
            grid_image.origin_col,
        ]
    )

    # Point to extrapolate in row col
    point = [0.5, 0.5]

    # Row array coordinates
    row_indexes = np.arange(grid_origin[0], grid_origin[0] + grid_image.nb_rows * grid_step[0], grid_step[0])
    # Col array coordinates
    col_indexes = np.arange(grid_origin[1], grid_origin[1] + grid_image.nb_columns * grid_step[1], grid_step[1])

    # Shareloc
    shareloc_outputs = gri_geotiff.direct_loc_h(point[0], point[1], 1000.0)

    # Extrapolation with scipy in 2D
    scipy_outputs = scipy.interpolate.interpn(
        (row_indexes, col_indexes), data_lon_1000, point, bounds_error=False, fill_value=None
    )

    # Test output equals to scipy ground truth
    np.testing.assert_allclose(shareloc_outputs[0][1], scipy_outputs[0])
