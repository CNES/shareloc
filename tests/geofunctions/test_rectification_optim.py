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
Test module for rectification grid interpolation class shareloc/geofunctions/rectification*.py
Ground truth references (otb_{left/right}_grid*.tif) have been generated using OTB StereoRectificationGridGenerator
application.
"""
# Standard imports
import math
import os

# Third party imports
import numpy as np
import pytest

# Shareloc bindings
import bindings_cpp
from shareloc.geofunctions.rectification import compute_local_epipolar_line, moving_along_axis

# Shareloc test imports
from ..helpers import DTMIntersection_constructor, bindings_cpp_constructor, data_path


@pytest.mark.unit_tests
def test_epipolar_angle():
    """
    test epipolar angle computation
    """
    # First case : same column, positive direction [row, col, alt]
    start_line_1 = [1, 0, 0]
    end_line_1 = [2, 0, 0]

    reference_alpha_1 = math.pi / 2.0
    alpha = bindings_cpp.compute_epipolar_angle(end_line_1, start_line_1)

    assert alpha == reference_alpha_1

    # Second case : same column, negative direction [row, col, alt]
    start_line_2 = [2, 0, 0]
    end_line_2 = [1, 0, 0]

    reference_alpha_2 = -(math.pi / 2.0)
    alpha = bindings_cpp.compute_epipolar_angle(end_line_2, start_line_2)

    assert alpha == reference_alpha_2

    # Third case : different column, positive direction [row, col, alt]
    start_line_3 = [2, 0, 0]
    end_line_3 = [1, 1, 0]

    slope = (1 - 2) / (1 - 0)

    reference_alpha_3 = np.arctan(slope)
    alpha = bindings_cpp.compute_epipolar_angle(end_line_3, start_line_3)

    assert alpha == reference_alpha_3

    # Fourth case : different column, negative direction [row, col, alt]
    start_line_4 = [2, 1, 0]
    end_line_4 = [1, 0, 0]

    slope = (1 - 2) / (0 - 1)
    reference_alpha_4 = math.pi + np.arctan(slope)
    alpha = bindings_cpp.compute_epipolar_angle(end_line_4, start_line_4)

    assert alpha == reference_alpha_4

    # With multiple point
    start_lines = np.stack((start_line_1, start_line_2, start_line_3, start_line_4))
    end_lines = np.stack((end_line_1, end_line_2, end_line_3, end_line_4))
    reference_alphas = np.stack((reference_alpha_1, reference_alpha_2, reference_alpha_3, reference_alpha_4))

    alphas = np.empty(4)
    for i in range(np.shape(start_lines)[0]):
        alphas[i] = bindings_cpp.compute_epipolar_angle(end_lines[i], start_lines[i])

    np.testing.assert_array_equal(alphas, reference_alphas)


# @pytest.mark.unit_tests
# def test_epipolar_angle():
#     """
#     test epipolar angle computation
#     """
#     # First case : same column, positive direction [row, col, alt]
#     start_line_1 = np.array([1, 0, 0])
#     end_line_1 = np.array([2, 0, 0])

#     reference_alpha_1 = math.pi / 2.0
#     alpha = rpc_c.compute_epipolar_angle(
#         [end_line_1[0]], [end_line_1[1]], [end_line_1[2]], [start_line_1[0]], [start_line_1[1]], [start_line_1[2]]
#     )

#     assert alpha[0] == reference_alpha_1

#     # Second case : same column, negative direction [row, col, alt]
#     start_line_2 = np.array([2, 0, 0])
#     end_line_2 = np.array([1, 0, 0])

#     reference_alpha_2 = -(math.pi / 2.0)
#     alpha = rpc_c.compute_epipolar_angle(
#         [end_line_2[0]], [end_line_2[1]], [end_line_2[2]], [start_line_2[0]], [start_line_2[1]], [start_line_2[2]]
#     )

#     assert alpha[0] == reference_alpha_2

#     # Third case : different column, positive direction [row, col, alt]
#     start_line_3 = np.array([2, 0, 0])
#     end_line_3 = np.array([1, 1, 0])

#     slope = (1 - 2) / (1 - 0)

#     reference_alpha_3 = np.arctan(slope)
#     alpha = rpc_c.compute_epipolar_angle(
#         [end_line_3[0]], [end_line_3[1]], [end_line_3[2]], [start_line_3[0]], [start_line_3[1]], [start_line_3[2]]
#     )

#     assert alpha == reference_alpha_3

#     # Fourth case : different column, negative direction [row, col, alt]
#     start_line_4 = np.array([2, 1, 0])
#     end_line_4 = np.array([1, 0, 0])

#     slope = (1 - 2) / (0 - 1)
#     reference_alpha_4 = math.pi + np.arctan(slope)
#     alpha = rpc_c.compute_epipolar_angle(
#         [end_line_4[0]], [end_line_4[1]], [end_line_4[2]], [start_line_4[0]], [start_line_4[1]], [start_line_4[2]]
#     )

#     assert alpha == reference_alpha_4

#     # With multiple point
#     start_lines = np.stack((start_line_1, start_line_2, start_line_3, start_line_4))
#     end_lines = np.stack((end_line_1, end_line_2, end_line_3, end_line_4))
#     reference_alphas = np.stack((reference_alpha_1, reference_alpha_2, reference_alpha_3, reference_alpha_4))

#     alphas = compute_epipolar_angle(end_lines, start_lines)
#     alphas = rpc_c.compute_epipolar_angle(
#         end_lines[:, 0], end_lines[:, 1], end_lines[:, 2], start_lines[:, 0], start_lines[:, 1], start_lines[:, 2]
#     )

#     np.testing.assert_array_equal(alphas, reference_alphas)


def test_rectification_moving_along_lines():
    """
    Test moving along line in epipolar geometry
    """

    geom_left_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "left_image.geom"))
    geom_right_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "right_image.geom"))

    current_coords = [5000.5, 5000.5, 0.0]
    mean_spacing = 1.0
    epi_step = 1
    alphas = 0.0
    default_elev = 0.0
    # ground truth next pixel
    # col pixel size of the image
    col_pixel_size = 1.0
    reference_next_cords = np.array([5000.5, 5000.5 + col_pixel_size, 0.0], dtype=np.float64)

    res_cpp = bindings_cpp.moving_along_axis(
        geom_left_cpp,
        geom_right_cpp,
        current_coords,
        mean_spacing,
        default_elev,
        epi_step,
        alphas,
        1,
    )

    assert res_cpp[0][0] == reference_next_cords[0]
    assert res_cpp[0][1] == reference_next_cords[1]
    assert res_cpp[0][2] == reference_next_cords[2]


def test_rectification_moving_to_next_lines():
    """
    Test moving to next line in epipolar geometry
    """

    geom_left_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "left_image.geom"))
    geom_right_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "right_image.geom"))

    current_coords = [5000.5, 5000.5, 0.0]
    mean_spacing = 1.0
    epi_step = 1
    alphas = 0.0
    default_elev = 0.0
    # ground truth next pixel
    # col pixel size of the image
    row_pixel_size = 1.0
    reference_next_cords = np.array([5000.5 + row_pixel_size, 5000.5, 0.0], dtype=np.float64)

    res_cpp = bindings_cpp.moving_along_axis(
        geom_left_cpp,
        geom_right_cpp,
        current_coords,
        mean_spacing,
        default_elev,
        epi_step,
        alphas,
        0,
    )

    assert res_cpp[0][0] == reference_next_cords[0]
    assert res_cpp[0][1] == reference_next_cords[1]
    assert res_cpp[0][2] == reference_next_cords[2]


def test_moving_along_axis(init_rpc_geom_model):
    """
    Test moving along axis function
    """

    geom_left_py, geom_right_py = init_rpc_geom_model
    geom_left_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "left_image.geom"))
    geom_right_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "right_image.geom"))

    mnt = os.path.join(data_path(), "dtm/srtm_ventoux/srtm90_non_void_filled/N44E005.hgt")

    dtm_py, dtm_cpp = DTMIntersection_constructor(mnt)

    current_coords = [5000.5, 5000.5, 0.0]
    current_coords_py = np.array([[5000.5, 5000.5, 0.0]], dtype=np.float64)
    mean_spacing = 1.0
    epi_step = 1
    alphas = 0.0

    res_cpp = bindings_cpp.moving_along_axis(
        geom_left_cpp,
        geom_right_cpp,
        current_coords,
        mean_spacing,
        dtm_cpp,
        epi_step,
        alphas,
        0,
    )

    res_left_py, res_rigth_py = moving_along_axis(
        geom_left_py,
        geom_right_py,
        current_coords_py,
        mean_spacing,
        dtm_py,
        epi_step,
        alphas,
        0,
    )

    assert res_cpp[0][0] == res_left_py[0][0]
    assert res_cpp[0][1] == res_left_py[0][1]
    assert res_cpp[0][2] == res_left_py[0][2]
    assert res_cpp[1][0] == res_rigth_py[0][0]
    assert res_cpp[1][1] == pytest.approx(res_rigth_py[0][1], abs=8e-12)
    assert res_cpp[1][2] == res_rigth_py[0][2]

    with pytest.raises(RuntimeError):
        res_cpp = bindings_cpp.moving_along_axis(
            geom_left_cpp,
            geom_right_cpp,
            current_coords,
            mean_spacing,
            dtm_cpp,
            epi_step,
            alphas,
            2,
        )


def test_compute_local_epipolar_line(init_rpc_geom_model):
    """
    Test compute_local_epipolar_line methode
    """

    geom_left_py, geom_right_py = init_rpc_geom_model
    geom_left_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "left_image.geom"))
    geom_right_cpp = bindings_cpp_constructor(os.path.join(data_path(), "rectification", "right_image.geom"))

    mnt = os.path.join(data_path(), "dtm/srtm_ventoux/srtm90_non_void_filled/N44E005.hgt")
    dtm_py, dtm_cpp = DTMIntersection_constructor(mnt)

    # with dtm
    left_point = np.array([140.3, 154.3, 50])
    res_py_0, res_py_1 = compute_local_epipolar_line(geom_left_py, geom_right_py, left_point, dtm_py, 50)
    res_cpp_0, res_cpp_1 = bindings_cpp.compute_local_epipolar_line(
        geom_left_cpp, geom_right_cpp, left_point, dtm_cpp, 50
    )
    np.testing.assert_allclose(res_py_0[0], res_cpp_0, 0, 2e-11)
    np.testing.assert_allclose(res_py_1[0], res_cpp_1, 0, 4e-12)

    # withm altitude
    left_point = np.array([140.3, 154.3, 50])
    res_py_0, res_py_1 = compute_local_epipolar_line(geom_left_py, geom_right_py, left_point, 75, 50)
    res_cpp_0, res_cpp_1 = bindings_cpp.compute_local_epipolar_line(geom_left_cpp, geom_right_cpp, left_point, 75, 50)
    np.testing.assert_allclose(res_py_0[0], res_cpp_0, 0, 8e-12)
    np.testing.assert_allclose(res_py_1[0], res_cpp_1, 0, 4e-12)
