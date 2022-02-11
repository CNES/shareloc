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
This module contains the rectification_grids class to handle CARS rectification grids.
"""
# pylint: disable=no-member

# Third party imports
import numpy as np
import rasterio as rio
from scipy import interpolate


class RectificationGrid:
    """rectification grid"""

    def __init__(self, grid_filename):
        """
        Constructor
        :param grid_filename: grid filename
        :type grid_filename: string
        """
        self.filename = grid_filename

        dataset = rio.open(grid_filename)

        transform = dataset.transform
        step_col = transform[0]
        step_row = transform[4]
        # 0 or 0.5
        [ori_col, ori_row] = transform * (0.5, 0.5)  # positions au centre pixel

        # print("ori {} {} step {} {}".format(ori_col,ori_y,step_x,step_y))
        last_col = ori_col + step_col * dataset.width
        last_row = ori_row + step_row * dataset.height

        # transform dep to positions
        self.row_dep = dataset.read(2).transpose()
        self.col_dep = dataset.read(1).transpose()

        cols = np.arange(ori_col, last_col, step_col)
        rows = np.arange(ori_row, last_row, step_row)
        self.grid_row, self.grid_col = np.mgrid[ori_col:last_col:step_col, ori_row:last_row:step_row]
        self.points = (cols, rows)
        self.row_positions = self.row_dep + self.grid_col
        self.col_positions = self.col_dep + self.grid_row

    def get_positions(self):
        """
        return grid positions
        :return grid positions
        :rtype  numpy array, numpy array
        """
        return self.row_positions, self.col_positions

    def interpolate(self, positions):
        """
        interpolate position
        :param positions : positions to interpolate : array  Nx2 [col,row]
        :type positions: np.array
        :return interpolated positions : array  Nx2 [col,row]
        :rtype  np.array
        """
        interp_row = interpolate.interpn(
            self.points, self.row_positions, positions, method="linear", bounds_error=False, fill_value=None
        )
        interp_col = interpolate.interpn(
            self.points, self.col_positions, positions, method="linear", bounds_error=False, fill_value=None
        )
        interp_pos = np.stack((interp_col, interp_row)).transpose()
        return interp_pos
