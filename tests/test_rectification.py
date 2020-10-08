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
import os
import pytest
import numpy as np

from shareloc.rectification.rectification_grid import rectification_grid



    

def test_rectification_grid_interpolation_regular():
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    rectif_grid = rectification_grid(grid_filename)
    lig = np.arange(15, 1000, 10)
    col = np.arange(15, 1000, 10)
    #value at position [15,15]
    value_lig = np.sum(rectif_grid.lig_dep[0:2,0:2]) / 4.0 + 15.0
    value_col = np.sum(rectif_grid.col_dep[0:2,0:2]) / 4.0 + 15.0
    coords = rectif_grid.interpolate_regular(lig,col)
    assert(value_col == pytest.approx(coords[0,0,1],abs=1e-4))
    assert(value_lig == pytest.approx(coords[0,0,0],abs=1e-4))


@pytest.mark.unit_tests
def test_rectification_grid_interpolation():
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    rectif_grid = rectification_grid(grid_filename)
    #value at position [15,15]

    bary_diag_lig = rectif_grid.lig_dep[0,0] + rectif_grid.lig_dep[1,1]
    bary_diag_lig = bary_diag_lig / 2.0 + 15.0
    bary_diag_col = rectif_grid.col_dep[0,0] + rectif_grid.col_dep[1,1]
    bary_diag_col = bary_diag_col / 2.0 + 15.0

    grid_y_out, grid_x_out = np.mgrid[15:60:15, 15:60:15]
    sensor_positions = np.stack((grid_x_out, grid_y_out), axis=2)
    coords = rectif_grid.interpolate(sensor_positions)
    assert(bary_diag_col == pytest.approx(coords[0,0,1],abs=1e-4))
    assert(bary_diag_lig == pytest.approx(coords[0,0,0],abs=1e-4))


