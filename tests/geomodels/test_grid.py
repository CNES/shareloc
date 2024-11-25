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
import os

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


def test_read_grid_ref():
    """
    The purpose of this test is to test grid handling with non EPSG based metadata REF key
    """
    grid_path = os.path.join(data_path(), "grid", "grid_ref_proj4.tif")
    grid = GeoModel(grid_path, "GRID")
    assert grid.epsg == 4326
    assert grid.repter == "+proj=latlon +R=1737400.0 +no_defs +type=crs"


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
def test_compare_extrapolation():
    """
    The purpose of this test is to give an idea of
    the extrapolation precision for PHR ventoux grid example ([200,200] line/pixels step) using two grids, which have
    different size
    """
    grid_path = os.path.join(data_path(), "grid", "phr_ventoux")
    gri = GeoModel(os.path.join(grid_path, "GRID_PHR1B_P_201308051042194_SEN_690908101-001.tif"), "GRID")
    gri_large = GeoModel(os.path.join(grid_path, "GRID_PHR1B_P_201308051042194_SEN_690908101-001_LARGE.tif"), "GRID")
    # Point to extrapolate in row col
    points = np.array(
        [[6000.5, 6000.5], [7100.5, 6000.5], [8000.5, 6000.5], [6000.5, 8100.5], [6000.5, 9000.5], [8000.5, 9000.5]]
    )
    lonlatalt = gri_large.direct_loc_h(points[:, 0], points[:, 1], 1000.0)
    gri_large.estimate_inverse_loc_predictor()
    row_inv, col_inv, _ = gri_large.inverse_loc(lonlatalt[:, 0], lonlatalt[:, 1], lonlatalt[:, 2])

    # 3 coordinates are inside the large grid, very close results are expected (1e-8 pixels)
    np.testing.assert_allclose(np.array((row_inv, col_inv)).transpose(), points, rtol=0.0, atol=1e-8)

    lonlatalt = gri_large.direct_loc_h(points[:, 0], points[:, 1], 1000.0)
    gri.estimate_inverse_loc_predictor()
    row_inv, col_inv, _ = gri.inverse_loc(lonlatalt[:, 0], lonlatalt[:, 1], lonlatalt[:, 2])
    inv_loc = np.array((row_inv, col_inv)).transpose()
    # first point is inside the grid
    np.testing.assert_allclose(inv_loc[0, :], points[0, :], rtol=0.0, atol=1e-8)

    # extrapolation of 100 rows
    np.testing.assert_allclose(inv_loc[1, :], points[1, :], rtol=0.0, atol=6e-3)

    # extrapolation of 1000 rows
    np.testing.assert_allclose(inv_loc[2, :], points[2, :], rtol=0.0, atol=6e-2)

    # extrapolation of 100 cols
    np.testing.assert_allclose(inv_loc[3, :], points[3, :], rtol=0.0, atol=8e-3)

    # extrapolation of 1000 cols
    np.testing.assert_allclose(inv_loc[4, :], points[4, :], rtol=0.0, atol=8e-2)

    # extrapolation of 1000 rows 1000 cols
    np.testing.assert_allclose(inv_loc[5, :], points[5, :], rtol=0.0, atol=4e-2)


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
