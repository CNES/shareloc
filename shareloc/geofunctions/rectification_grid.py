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
    """
    Rectification grid
    """

    def __init__(self, grid_filename: str, is_displacement_grid: bool = False, interpolator: str = "linear"):
        """
        Constructor
        :param grid_filename: grid filename
        :type filename: string
        :param is_displacement_grid: True if is a displacement grid
        :type is_displacement_grid: bool
        :param interpolator: grid interpolator for scipy/interpolate.RegularGridInterpolator
        :type interpolator: str
        """
        self.grid_filename = grid_filename

        dataset = rio.open(grid_filename)

        transform = dataset.transform
        step_col = transform[0]
        step_row = transform[4]
        # 0 or 0.5
        [ori_col, ori_row] = transform * (0.5, 0.5)  # center pixel position

        # print("ori {} {} step {} {}".format(ori_col,ori_y,step_x,step_y))
        last_col = ori_col + step_col * dataset.width
        last_row = ori_row + step_row * dataset.height

        cols = np.arange(ori_col, last_col, step_col)
        rows = np.arange(ori_row, last_row, step_row)
        self.points = (cols, rows)

        if is_displacement_grid:
            # transform dep to positions
            self.row_dep = dataset.read(2).transpose()
            self.col_dep = dataset.read(1).transpose()

            self.grid_row, self.grid_col = np.mgrid[ori_col:last_col:step_col, ori_row:last_row:step_row]
            self.points = (cols, rows)
            self.row_positions = self.row_dep + self.grid_col
            self.col_positions = self.col_dep + self.grid_row
        else:
            self.row_positions = dataset.read(2).transpose()
            self.col_positions = dataset.read(1).transpose()

        self.interpolator = interpolate.RegularGridInterpolator(
            self.points,
            np.stack((self.col_positions, self.row_positions), axis=2),
            method=interpolator,
            bounds_error=False,
            fill_value=None,
        )

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

        return self.interpolator(positions)
