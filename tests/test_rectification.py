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



    
@pytest.mark.parametrize("row,col", [(15,0)])
@pytest.mark.unit_tests
def test_rectification_grid_interpolation_one_point(row,col):
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    rectif_grid = rectification_grid(grid_filename)
    #value at position [15,15]
    value_row = np.sum(rectif_grid.row_positions[0,0:2]) / 2.0
    value_col = np.sum(rectif_grid.col_positions[0,0:2]) / 2.0
    coords = rectif_grid.interpolate((col,row))
    assert(value_row == pytest.approx(coords[0,1],abs=1e-4))
    assert(value_col == pytest.approx(coords[0,0],abs=1e-4))


@pytest.mark.unit_tests
def test_rectification_grid_interpolation():
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    rectif_grid = rectification_grid(grid_filename)
    #value at position [15,15]

    value_row = np.sum(rectif_grid.row_positions[0:2,0:2]) / 4.0
    value_col = np.sum(rectif_grid.col_positions[0:2,0:2]) / 4.0
    sensor_positions = np.zeros((2,2))
    sensor_positions[0,:] = [15.0,15.0]
    sensor_positions[1,:] = [0.0,0.0]
    coords = rectif_grid.interpolate(sensor_positions)
    assert(value_col == pytest.approx(coords[0,0],abs=1e-4))
    assert(value_row == pytest.approx(coords[0,1],abs=1e-4))


@pytest.mark.unit_tests
def test_rectification_grid_extrapolation():
    """
    Test interpolation on rectification grid
    """
    grid_filename = os.path.join(os.environ["TESTPATH"], "rectification_grids", "left_epipolar_grid_ventoux.tif")
    rectif_grid = rectification_grid(grid_filename)

    sensor_positions = np.zeros((2,2))
    sensor_positions[0,:] = [30.0,-10.0]
    sensor_positions[1,:] = [30.0,691.0]
    coords = rectif_grid.interpolate(sensor_positions)
    assert(5065.72347005208303016 == pytest.approx(coords[0,1],abs=1e-10))
    assert(4883.84894205729142413474619389 == pytest.approx(coords[1,1],abs=1e-10))
