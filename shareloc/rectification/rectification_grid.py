
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

        lig_dep = dataset.read(2)
        col_dep = dataset.read(1)
        epipolar_shape =  lig_dep.shape
        x = np.arange(ori_x, last_x, step_x)
        y = np.arange(ori_y, last_y, step_y)

        #transform dep to positions
        grid_x_out, grid_y_out = np.mgrid[ori_y:last_y:step_y, ori_x:last_x:step_x]
        self.lig_positions = lig_dep + grid_x_out
        self.col_positions = col_dep + grid_y_out

        self.interp_lig = interpolate.interp2d(x, y, self.lig_positions, kind='linear')
        self.interp_col = interpolate.interp2d(x, y, self.col_positions, kind='linear')


    def interpolate_regular(self,lig,col):
        interp_lig = self.interp_lig(lig,col)
        interp_col = self.interp_col(lig,col)
        coord = np.stack((interp_lig,interp_col), axis = 2)
        return coord

    def interpolate_point_list(self,point_list):
        return None

    #gridata