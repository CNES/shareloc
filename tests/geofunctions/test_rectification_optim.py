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

# Third party imports
import numpy as np
import pytest

import bindings_cpp


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
