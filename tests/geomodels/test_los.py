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
Module to test functions that use los
"""

import os

import numpy as np
import pytest

# Shareloc imports
from shareloc.geomodels import GeoModel
from shareloc.geomodels.los import LOS
from shareloc.proj_utils import coordinates_conversion

# Shareloc test imports
from ..helpers import data_path


@pytest.mark.parametrize("alt", [None, ([310, 850])])
def test_los_creation_with_different_alt(alt):
    """
    Test los_creation
    """

    data_folder = data_path()

    # Load matches
    matches = np.load(os.path.join(data_folder, "triangulation/matches-crop.npy"))
    matches_left = matches[:, 0:2]

    # Load geometrical model
    id_scene = "P1BP--2017092838284574CP"
    file_dimap = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    geometrical_model = GeoModel(file_dimap)
    # Get alt max
    if alt is None:
        alt_min, alt_max = geometrical_model.get_alt_min_max()
    else:
        alt_min, alt_max = alt

    # Create los
    left_los = LOS(matches_left, geometrical_model, alt)

    # Create ground truth
    los_nb_gt = 121

    # LOS construction
    los_extrema = np.zeros([2 * los_nb_gt, 3])
    list_col, list_row = (matches_left[:, 0], matches_left[:, 1])
    los_extrema[np.arange(0, 2 * los_nb_gt, 2), :] = geometrical_model.direct_loc_h(list_row, list_col, alt_max)
    los_extrema[np.arange(1, 2 * los_nb_gt, 2), :] = geometrical_model.direct_loc_h(list_row, list_col, alt_min)

    in_crs = 4326
    out_crs = 4978
    ecef_coord = coordinates_conversion(los_extrema, in_crs, out_crs)
    sis_gt = ecef_coord[0::2, :]

    vis_gt = sis_gt - ecef_coord[1::2, :]
    # /!\ normalized
    #
    #  direction vector creation
    vis_norm = np.linalg.norm(vis_gt, axis=1)
    rep_vis_norm = np.tile(vis_norm, (3, 1)).transpose()
    vis_gt = vis_gt / rep_vis_norm

    assert left_los.number == los_nb_gt
    assert (left_los.starting_points == sis_gt).all
    assert (left_los.viewing_vectors == vis_gt).all


def test_compare_two_altitude():
    """
    Test los_creation with two different altitudes
    and test that sis and vis aren't equal
    """

    data_folder = data_path()

    # Load matches
    matches = np.load(os.path.join(data_folder, "triangulation/matches-crop.npy"))
    matches_left = matches[:, 0:2]

    # Load geometrical model
    id_scene = "P1BP--2017092838284574CP"
    file_dimap = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    geometrical_model = GeoModel(file_dimap)

    # Create los
    model_los = LOS(matches_left, geometrical_model, [310, 850])
    wrong_los = LOS(matches_left, geometrical_model, [0, 210])

    assert (model_los.starting_points != wrong_los.starting_points).all
    assert (model_los.viewing_vectors != wrong_los.viewing_vectors).all


def test_los_parameters():
    """
    test los parameters
    """

    data_folder = data_path()

    # Load matches
    matches = np.load(os.path.join(data_folder, "triangulation/matches-crop.npy"))
    matches_left = matches[:, 0:2]

    # Load geometrical model
    id_scene = "P1BP--2017092838284574CP"
    file_dimap = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    geometrical_model = GeoModel(file_dimap)

    # Create los
    model_los = LOS(matches_left, geometrical_model, [310, 850])

    np.testing.assert_allclose(
        model_los.starting_points[0],
        np.array([4581872.886817954, 566896.0588669669, 4387122.994634065]),
        rtol=0,
        atol=0,
    )
    np.testing.assert_allclose(
        model_los.viewing_vectors[0],
        np.array([0.6168861918045335, 0.0018007411543079202, 0.7870503056934771]),
        rtol=0,
        atol=0,
    )
    np.testing.assert_allclose(
        model_los.ending_points[0], np.array([4581535.247284677, 566895.0732696055, 4386692.2193931695]), rtol=0, atol=0
    )
    assert model_los.number == 121
