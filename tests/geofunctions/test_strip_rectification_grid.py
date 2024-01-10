#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).
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
Test module for rectification grid interpolation class shareloc/geofunctions/rectification*.py
"""

# Standard imports
import os

import numpy as np

# Third party imports
import pytest
import rasterio

# Shareloc imports
from shareloc.geofunctions.rectification import compute_strip_of_epipolar_grid, positions_to_displacement_grid

# Shareloc test imports
from tests.helpers import data_path


@pytest.mark.unit_tests
def test_compute_strip_of_epipolar_grid_columns_lines_rectangular(
    init_inputs_rectification_fixture, init_rpc_geom_model
):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio
    test with non squared grid

    Input Geomodels: RPC
    Earth elevation: default to 0.0
    """
    geom_model_left, geom_model_right = init_rpc_geom_model

    (
        left_position_point,
        right_position_point,
        spacing,
        epi_step,
        elevation_offset,
        default_elev,
        grid_size,
    ) = init_inputs_rectification_fixture

    # Change size to make the grid rectangular
    grid_size[0] = 21

    left_grid, right_grid, alphas, mean_br_col = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_position_point,
        right_position_point,
        spacing,
        axis=0,
        strip_size=grid_size[0],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
    )

    left_grid, right_grid, alphas, mean_br = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_grid,
        right_grid,
        spacing,
        axis=1,
        strip_size=grid_size[1],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
        epipolar_angles=alphas,
    )
    mean_br = (mean_br * (grid_size[1] * (grid_size[0] - 1)) + mean_br_col * grid_size[0]) / (
        grid_size[1] * grid_size[0]
    )

    # gt reference
    reference_left_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_left.tif")
    ).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_right.tif")
    ).read()

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1][:21, :], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0][:21, :], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1][:21, :], right_grid[:, :, 0])
    np.testing.assert_allclose(reference_right_grid[0][:21, :], right_grid[:, :, 1], atol=2.0e-12)

    # Check mean_baseline_ratio
    reference_mean_br = 0.7024809
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)


@pytest.mark.unit_tests
def test_compute_strip_of_epipolar_grid_columns_lines(init_inputs_rectification_fixture, init_rpc_geom_model):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: default to 0.0
    """
    geom_model_left, geom_model_right = init_rpc_geom_model

    (
        left_position_point,
        right_position_point,
        spacing,
        epi_step,
        elevation_offset,
        default_elev,
        grid_size,
    ) = init_inputs_rectification_fixture

    left_grid, right_grid, alphas, mean_br_col = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_position_point,
        right_position_point,
        spacing,
        axis=0,
        strip_size=grid_size[0],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
    )

    left_grid, right_grid, alphas, mean_br = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_grid,
        right_grid,
        spacing,
        axis=1,
        strip_size=grid_size[1],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
        epipolar_angles=alphas,
    )
    mean_br = (mean_br * (grid_size[1] * (grid_size[0] - 1)) + mean_br_col * grid_size[0]) / (
        grid_size[1] * grid_size[0]
    )

    # shareloc reference
    reference_left_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_left.tif")
    ).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_right.tif")
    ).read()

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1], right_grid[:, :, 0])
    np.testing.assert_allclose(reference_right_grid[0], right_grid[:, :, 1], atol=2.0e-12)

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio
    reference_mean_br = 0.7040047235162911
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-6)


@pytest.mark.unit_tests
def test_compute_strip_of_epipolar_grid_lines_columns(init_inputs_rectification_fixture, init_rpc_geom_model):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: default to 0.0
    """
    geom_model_left, geom_model_right = init_rpc_geom_model

    (
        left_position_point,
        right_position_point,
        spacing,
        epi_step,
        elevation_offset,
        default_elev,
        grid_size,
    ) = init_inputs_rectification_fixture

    left_grid, right_grid, alphas, mean_br_line = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_position_point,
        right_position_point,
        spacing,
        axis=1,
        strip_size=grid_size[1],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
    )

    left_grid, right_grid, _, mean_br = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_grid,
        right_grid,
        spacing,
        axis=0,
        strip_size=grid_size[0],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
        epipolar_angles=alphas,
    )
    mean_br = (mean_br * (grid_size[1] * (grid_size[0] - 1)) + mean_br_line * grid_size[0]) / (
        grid_size[1] * grid_size[0]
    )

    # gt reference
    reference_left_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_left.tif")
    ).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_right.tif")
    ).read()

    # Check epipolar grids
    assert reference_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=7e-3)
    assert reference_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=7e-3)

    assert reference_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=7e-3)
    assert reference_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=7e-3)

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio
    reference_mean_br = 0.7040047235162911
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)


@pytest.mark.unit_tests
def test_positions_to_displacement_grid():
    """
    Test displacement grids generation : check epipolar grids

    Input position grids
    """

    position_grid_left = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_left.tif")
    ).read()
    position_grid_right = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_right.tif")
    ).read()
    position_grid_left = np.moveaxis(position_grid_left, 0, 2)[:, :, [1, 0]]
    position_grid_right = np.moveaxis(position_grid_right, 0, 2)[:, :, [1, 0]]

    left_grid, right_grid, grid_geotransform = positions_to_displacement_grid(
        position_grid_left, position_grid_right, epi_step=30.0
    )

    # Shareloc reference
    reference_left_dataset = rasterio.open(os.path.join(data_path(), "rectification", "shareloc_gt_left_grid.tif"))
    reference_left_grid = reference_left_dataset.read()
    reference_right_dataset = rasterio.open(os.path.join(data_path(), "rectification", "shareloc_gt_right_grid.tif"))
    reference_right_grid = reference_right_dataset.read()

    assert grid_geotransform == reference_left_dataset.transform
    assert grid_geotransform == reference_right_dataset.transform

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1], right_grid[:, :, 0])
    np.testing.assert_array_equal(reference_right_grid[0], right_grid[:, :, 1])


# pylint: disable=too-many-locals
def test_working_with_strip(init_rpc_geom_model, init_inputs_rectification_fixture):
    """
    Test epipolar creation grid by strip. This test presents how to use computation by strips.
    It is useful for describing parralelisation opportunity
    """
    geom_model_left, geom_model_right = init_rpc_geom_model
    (
        left_position_point,
        right_position_point,
        spacing,
        epi_step,
        elevation_offset,
        default_elev,
        grid_size,
    ) = init_inputs_rectification_fixture

    # choose by user
    axis = 0
    size_strip = 8

    # First row is computed
    first_left_points, first_right_points, first_alphas, mean_br_col = compute_strip_of_epipolar_grid(
        geom_model_left,
        geom_model_right,
        left_position_point,
        right_position_point,
        spacing,
        axis=0,
        strip_size=grid_size[0],
        epi_step=epi_step,
        elevation=default_elev,
        elevation_offset=elevation_offset,
    )
    sum_br = mean_br_col * grid_size[0]

    left_strips = []
    right_strips = []
    alphas_strip = []

    # create strips
    for i in range(0, grid_size[axis], size_strip):
        # work only if first used axis is 0
        # else strip = first_points[:, :, i:i+size_strip_columns, :]
        left_strip = first_left_points[i : i + size_strip, :, :]
        right_strip = first_right_points[i : i + size_strip, ::]
        alphas = first_alphas[i : i + size_strip, :]
        left_strips.append(left_strip)
        right_strips.append(right_strip)
        alphas_strip.append(alphas)

    # instantiate second axis
    axis = 1

    # instantiate outputs
    left_grid, right_grid = np.zeros((grid_size[0], grid_size[1], 2)), np.zeros((grid_size[0], grid_size[1], 2))

    # loop on the strip
    for idx, left_strip in enumerate(left_strips):
        strip_left_grid, strip_right_grid, _, mbr = compute_strip_of_epipolar_grid(
            geom_model_left,
            geom_model_right,
            left_strip,
            right_strips[idx],
            spacing,
            axis=1,
            strip_size=grid_size[1],
            epi_step=epi_step,
            elevation=default_elev,
            elevation_offset=elevation_offset,
            epipolar_angles=alphas_strip[idx],
        )

        left_grid[idx * size_strip : size_strip + idx * size_strip, :, :] = strip_left_grid[:, :, 0:2]
        right_grid[idx * size_strip : size_strip + idx * size_strip, :, :] = strip_right_grid[:, :, 0:2]

        sum_br = sum_br + mbr * left_strip.shape[0] * (grid_size[1] - 1)

    # Shareloc reference
    reference_left_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_left.tif")
    ).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_positions_grid_right.tif")
    ).read()

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1], right_grid[:, :, 0])
    np.testing.assert_allclose(reference_right_grid[0], right_grid[:, :, 1], atol=2.0e-12)

    mean_br = sum_br / (grid_size[0] * grid_size[1])
    reference_mean_br = 0.704004723
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-7)
