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

#-------------------------------------------------------------------------------

import numpy as np


class Sensor:
    """ base class for localization function.
    Underlying model can be both multi layer localization grids or RPCs models
    """
    def __init__(self, grid = None, rpc = None, dtm = None):
        """
        If grid and RPC are set grid model is used
        :param shareloc.grid grid: multi layer grid object
        :param shareloc.rpc rpc: rpc object
        :param shareloc.dtm dtm: dtm model
        """
        self.dtm = dtm
        self.grid = grid
        self.rpc = rpc
        self.use_rpc = False
        if self.rpc is not None and self.grid is None:
            self.use_rpc = True


    def forward(self, row, col, h):
        """
        forward localization
        :param row :  sensor row
        :param col : sensor col
        :param h: altitude
        :return coordinates : [lon,lat,h] (3D np.array)
        """
        if self.use_rpc == True:
            print('use rpc')
            coord = np.zeros(3)
            coord[0:2] = self.rpc.evalue_loc_d(col,row, h)
            coord[2] = h
            return coord
        else:
            return self.grid.fct_locdir_h(row,col,h)

    def forward_dtm(self, row, col):
        """
        forward localization on dtm
        :param row : sensor row
        :param col : sensor col
        :return coordinates : [lon,lat,h] (3D np.array)
        """
        if self.use_rpc == True:
            print('forward_dtm not yet impelemented for RPC model')
            return None
        else:
            if self.dtm is not None:
                return self.grid.fct_locdir_dtm(row,col, self.dtm)
            else:
                print('forward_dtm needs a dtm')
                return None

    def inverse(self,lon,lat,h):
        """
        inverse localization
        :param lat :  latitude
        :param lon : longitude
        :param h : altitude
        :return coordinates : [row,col,valid] (2D np.array), valid == 1 if coordinates is valid
        :rtype numpy.array
        """
        if self.use_rpc == True:
            print('use rpc')
            [col, row] = self.rpc.evalue_loc_i(lon,lat,h)
            return row,col,1
        else:
            if not hasattr(self.grid, 'pred_ofset_scale_lon'):
                self.grid.init_pred_loc_inv()
            return self.grid.fct_locinv([lon,lat,h])



