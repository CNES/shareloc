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
from shareloc.geomodels import GeoModel
from shareloc.image import Image

# Shareloc test imports
from tests.helpers import data_path


@pytest.fixture(name="get_geotiff_grid")
def fixture_get_geotiff_grid():
    """
    get grid and associated path
    """
    geotiff_grid_path = data_path("ellipsoide", "loc_direct_grid_PHR_2013072139303958CP.tif")
    gri_geotiff = GeoModel(geotiff_grid_path, "GRID")
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
    point = np.array([[0.5, 0.5], [-0.5, -0.5], [490.0, 500.0]])

    # Row array coordinates
    row_indexes = np.arange(grid_origin[0], grid_origin[0] + grid_image.nb_rows * grid_step[0], grid_step[0])
    # Col array coordinates
    col_indexes = np.arange(grid_origin[1], grid_origin[1] + grid_image.nb_columns * grid_step[1], grid_step[1])

    # Shareloc
    shareloc_outputs = gri_geotiff.direct_loc_h(point[:, 0], point[:, 1], 1000.0)

    # Extrapolation with scipy in 2D
    scipy_outputs = scipy.interpolate.interpn(
        (row_indexes, col_indexes), data_lon_1000, point, bounds_error=False, fill_value=None
    )

    # Test output equals to scipy ground truth
    np.testing.assert_allclose(shareloc_outputs[:, 1], scipy_outputs, rtol=0.0, atol=1e-9)

    # inverse loc extrapolation
    gri_geotiff.estimate_inverse_loc_predictor()
    row_inv, col_inv, _ = gri_geotiff.inverse_loc(
        shareloc_outputs[:, 0], shareloc_outputs[:, 1], shareloc_outputs[:, 2]
    )
    np.testing.assert_allclose(np.array((row_inv, col_inv)).transpose(), point, rtol=0.0, atol=1e-7)


@pytest.mark.unit_tests
def test_sensor_loc_nan(get_geotiff_grid):
    """
    Test direct localization  with grid containing NaN
    """
    h = 0
    row = np.array([np.nan])
    col = np.array([np.nan])
    gri, _ = get_geotiff_grid
    loc = gri.direct_loc_h(row, col, h)
    assert np.array_equal(loc, np.array([[np.nan, np.nan, h]]), equal_nan=True)
    gri.estimate_inverse_loc_predictor()
    lon = np.array([np.nan])
    lat = np.array([np.nan])
    res_inv = gri.inverse_loc(lon, lat, h)
    np.testing.assert_allclose(res_inv, [[np.nan], [np.nan], [0]], rtol=0, atol=1e-9)


@pytest.mark.unit_tests
def test_sensor_loc_inv(get_geotiff_grid):
    """
    test grid readers
    """
    gri, _ = get_geotiff_grid
    lon = 2.1828713504608683
    lat = 48.942429997483146
    h = 200.0
    gri.estimate_inverse_loc_predictor()
    res = gri.inverse_loc(lon, lat, h)
    np.testing.assert_allclose(res, [[50], [100], [200.0]], rtol=0, atol=1e-9)
    lon = np.array([2.1828713504608683, np.nan, 2.1828713504608683, np.nan, 2.1828713504608683])
    lat = np.array([48.942429997483146, np.nan, np.nan, 48.942429997483146, 48.942429997483146])
    h = 200.0
    res = gri.inverse_loc(lon, lat, h)
    np.testing.assert_allclose(
        res,
        [[50, np.nan, np.nan, np.nan, 50], [100, np.nan, np.nan, np.nan, 100], [200.0, 200.0, 200.0, 200.0, 200.0]],
        rtol=0,
        atol=1e-9,
    )
