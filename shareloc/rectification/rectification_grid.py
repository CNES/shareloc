
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

import rasterio as rio
import numpy as np
from scipy import interpolate

class rectification_grid:
    """ rectification grid
    """
    def __init__(self, grid_filename):
        """
        Constructor
        :param grid_filename: grid filename
        :type grid_filename: string
        """
        self.filename = grid_filename

        dataset = rio.open(grid_filename)

        transform = dataset.transform
        step_x = transform[0]
        step_y = transform[4]
        # 0 or 0.5
        [ori_x, ori_y] = transform*(0.5, 0.5) # positions au centre pixel

        #print("ori {} {} step {} {}".format(ori_x,ori_y,step_x,step_y))
        last_x = ori_x + step_x * dataset.width
        last_y = ori_y + step_y * dataset.height

        #transform dep to positions
        self.row_dep = dataset.read(2).transpose()
        self.col_dep = dataset.read(1).transpose()

        #TODO change notations
        x = np.arange(ori_x, last_x, step_x)
        y = np.arange(ori_y, last_y, step_y)
        self.grid_y, self.grid_x = np.mgrid[ori_x:last_x:step_x, ori_y:last_y:step_y]
        self.points = (x,y)
        self.row_positions = self.row_dep + self.grid_x
        self.col_positions = self.col_dep + self.grid_y


    def interpolate(self, positions):
        """
        interpolate position
        :param positions : positions to interpolate : array  Nx2 [col,row]
        :type positions: np.array
        :return interpolated positions : array  Nx2 [col,row]
        :rtype  np.array
        """
        #interp with grid data
        interp_row = interpolate.interpn(self.points, self.row_positions, positions, method='linear')
        interp_col = interpolate.interpn(self.points, self.col_positions, positions, method='linear')
        interp_pos = np.stack((interp_col,interp_row)).transpose()
        return interp_pos


    #gridata