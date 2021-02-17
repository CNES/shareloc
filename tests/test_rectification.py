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
Test module for rectification grid interpolation class shareloc/rectification/rectification_grid.py
"""

import os
import math
import pytest
import numpy as np
import rasterio

from shareloc.rpc.rpc import RPC
from shareloc.rectification.rectification_grid import RectificationGrid
from shareloc.rectification.rectification import (
    compute_epipolar_angle,
    moving_along_lines,
    moving_to_next_line,
    prepare_rectification,
    compute_stereorectification_epipolar_grids,
)
from shareloc.image.image import Image


@pytest.mark.parametrize("row,col", [(15, 0)])
@pytest.mark.unit_tests
def test_rectification_grid_interpolation_one_point(row, col):
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(os.environ["TESTPATH"], "rectification_grids", "grid_{}.tif".format(id_scene_right))
    rectif_grid = RectificationGrid(grid_filename)
    # value at position [15,15]
    value_row = np.sum(rectif_grid.row_positions[0, 0:2]) / 2.0
    value_col = np.sum(rectif_grid.col_positions[0, 0:2]) / 2.0
    coords = rectif_grid.interpolate((col, row))
    assert value_row == pytest.approx(coords[0, 1], abs=1e-4)
    assert value_col == pytest.approx(coords[0, 0], abs=1e-4)


@pytest.mark.unit_tests
def test_rectification_grid_interpolation():
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(os.environ["TESTPATH"], "rectification_grids", "grid_{}.tif".format(id_scene_right))

    rectif_grid = RectificationGrid(grid_filename)
    # value at position [15,15]

    value_row = np.sum(rectif_grid.row_positions[0:2, 0:2]) / 4.0
    value_col = np.sum(rectif_grid.col_positions[0:2, 0:2]) / 4.0
    sensor_positions = np.zeros((2, 2))
    sensor_positions[0, :] = [15.0, 15.0]
    sensor_positions[1, :] = [0.0, 0.0]
    coords = rectif_grid.interpolate(sensor_positions)
    assert value_col == pytest.approx(coords[0, 0], abs=1e-4)
    assert value_row == pytest.approx(coords[0, 1], abs=1e-4)


@pytest.mark.unit_tests
def test_rectification_grid_extrapolation():
    """
    Test interpolation on rectification grid
    """
    grid_filename = os.path.join(os.environ["TESTPATH"], "rectification_grids", "left_epipolar_grid_ventoux.tif")
    rectif_grid = RectificationGrid(grid_filename)

    sensor_positions = np.zeros((2, 2))
    sensor_positions[0, :] = [30.0, -10.0]
    sensor_positions[1, :] = [30.0, 691.0]
    coords = rectif_grid.interpolate(sensor_positions)
    assert pytest.approx(coords[0, 1], abs=1e-10) == 5065.72347005208303016
    assert pytest.approx(coords[1, 1], abs=1e-10) == 4883.84894205729142413474619389


@pytest.mark.unit_tests
def test_prepare_rectification():
    """
    Test prepare rectification : check grids size, epipolar image size, and left epipolar starting point
    """
    left_im = Image(os.path.join(os.environ["TESTPATH"], "rectification", "left_image.tif"))

    geom_model_left = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "left_image.geom"), topleftconvention=True
    )
    geom_model_right = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "right_image.geom"), topleftconvention=True
    )

    epi_step = 30
    elevation_offset = 50
    __, grid_size, rectified_image_size, left_epi_origin = prepare_rectification(
        left_im, geom_model_left, geom_model_right, epi_step, elevation_offset
    )
    # check size of the epipolar grids
    assert grid_size[0] == 22
    assert grid_size[1] == 22

    #  check size of rectified images
    assert rectified_image_size[0] == 612
    assert rectified_image_size[1] == 612

    # check the first epipolar point in the left image
    # ground truth values from OTB
    otb_output_origin_in_left_image = [5625.78139593690139008685946465, 5034.15635707952696975553408265, 0]

    # OTB convention is [col, row, altitude], shareloc convention is [row, col, altitude]
    assert otb_output_origin_in_left_image[1] == pytest.approx(left_epi_origin[0], abs=1e-5)
    assert otb_output_origin_in_left_image[0] == pytest.approx(left_epi_origin[1], abs=1e-5)

    assert left_epi_origin[2] == otb_output_origin_in_left_image[2]


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids():
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio
    """
    left_im = Image(os.path.join(os.environ["TESTPATH"], "rectification", "left_image.tif"))
    right_im = Image(os.path.join(os.environ["TESTPATH"], "rectification", "right_image.tif"))

    geom_model_left = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "left_image.geom"), topleftconvention=True
    )
    geom_model_right = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "right_image.geom"), topleftconvention=True
    )

    epi_step = 30
    elevation_offset = 50
    dem = None
    default_elev = 0.0
    left_grid, right_grid, img_size_row, img_size_col, mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, None, default_elev, epi_step, elevation_offset
    )

    gt_left_grid = rasterio.open(os.path.join(os.environ["TESTPATH"], "rectification", "gt_left_grid.tif")).read()
    gt_right_grid = rasterio.open(os.path.join(os.environ["TESTPATH"], "rectification", "gt_right_grid.tif")).read()

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert gt_left_grid[1] == pytest.approx(left_grid.data[0, :, :], abs=1e-2)
    assert gt_left_grid[0] == pytest.approx(left_grid.data[1, :, :], abs=1e-2)

    assert gt_right_grid[1] == pytest.approx(right_grid.data[0, :, :], abs=1e-2)
    assert gt_right_grid[0] == pytest.approx(right_grid.data[1, :, :], abs=1e-2)

    #  Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    #  Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    gt_mean_br = 0.704004705
    assert mean_br == pytest.approx(gt_mean_br, abs=1e-5)


@pytest.mark.unit_tests
def test_rectification_moving_along_line():
    """
    Test moving along line in epipolar geometry
    """
    geom_model_left = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "left_image.geom"), topleftconvention=True
    )
    geom_model_right = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "right_image.geom"), topleftconvention=True
    )

    current_left_coords = np.array([[5000.5, 5000.5, 0.0]], dtype=np.float64)
    mean_spacing = 1
    epi_step = 1
    alphas = 0

    # ground truth next pixel
    # col pixel size of the image
    col_pixel_size = 1.0
    gt_next_cords = np.array([[5000.5, 5000.5 + col_pixel_size, 0.0]], dtype=np.float64)

    next_cords, _ = moving_along_lines(
        geom_model_left, geom_model_right, current_left_coords, mean_spacing, epi_step, alphas
    )

    np.testing.assert_array_equal(gt_next_cords, next_cords)


@pytest.mark.unit_tests
def test_rectification_moving_to_next_line():
    """
    Test moving to next line in epipolar geometry
    """
    geom_model_left = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "left_image.geom"), topleftconvention=True
    )
    geom_model_right = RPC.from_any(
        os.path.join(os.environ["TESTPATH"], "rectification", "right_image.geom"), topleftconvention=True
    )

    current_left_coords = np.array([5000.5, 5000.5, 0.0], dtype=np.float64)
    mean_spacing = 1
    epi_step = 1
    alphas = 0

    # ground truth next pixel
    # row pixel size of the image
    row_pixel_size = 1.0
    gt_next_cords = np.array([5000.5 + row_pixel_size, 5000.5, 0.0], dtype=np.float64)

    next_cords, _ = moving_to_next_line(
        geom_model_left, geom_model_right, current_left_coords, mean_spacing, epi_step, alphas
    )

    np.testing.assert_array_equal(gt_next_cords, next_cords)


@pytest.mark.unit_tests
def test_epipolar_angle():
    """
    test epipolar angle computation
    """
    # First case : same column, positive direction [row, col, alt]
    start_line_1 = np.array([1, 0, 0])
    end_line_1 = np.array([2, 0, 0])

    gt_alpha_1 = math.pi / 2.0
    alpha = compute_epipolar_angle(end_line_1, start_line_1)
    assert alpha == gt_alpha_1

    # Second case : same column, negative direction [row, col, alt]
    start_line_2 = np.array([2, 0, 0])
    end_line_2 = np.array([1, 0, 0])

    gt_alpha_2 = -(math.pi / 2.0)
    alpha = compute_epipolar_angle(end_line_2, start_line_2)
    assert alpha == gt_alpha_2

    # Third case : different column, positive direction [row, col, alt]
    start_line_3 = np.array([2, 0, 0])
    end_line_3 = np.array([1, 1, 0])

    slope = (1 - 2) / (1 - 0)

    gt_alpha_3 = np.arctan(slope)
    alpha = compute_epipolar_angle(end_line_3, start_line_3)
    assert alpha == gt_alpha_3

    # Fourth case : different column, negative direction [row, col, alt]
    start_line_4 = np.array([2, 1, 0])
    end_line_4 = np.array([1, 0, 0])

    slope = (1 - 2) / (0 - 1)
    gt_alpha_4 = math.pi + np.arctan(slope)
    alpha = compute_epipolar_angle(end_line_4, start_line_4)
    assert alpha == gt_alpha_4

    # With multiple point
    start_lines = np.stack((start_line_1, start_line_2, start_line_3, start_line_4))
    end_lines = np.stack((end_line_1, end_line_2, end_line_3, end_line_4))
    gt_alphas = np.stack((gt_alpha_1, gt_alpha_2, gt_alpha_3, gt_alpha_4))

    alphas = compute_epipolar_angle(end_lines, start_lines)
    np.testing.assert_array_equal(alphas, gt_alphas)
