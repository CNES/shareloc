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


class Localization:
    """ base class for localization function.
    Underlying model can be both multi layer localization grids or RPCs models
    """
    def __init__(self, model, dtm = None):
        """
        constructor
        :param model : geometric model
        :type model  : shareloc.grid or  shareloc.rpc
        :param dtm  : dtm (optional)
        :type dtm  : shareloc.dtm
        """
        self.dtm = dtm
        self.use_rpc = model.type is 'rpc'
        self.model = model


    def direct(self, row, col, h = None):
        """
        direct localization
        :param row :  sensor row
        :type row : float
        :param col : sensor col
        :type col : float
        :param h: altitude, if none DTM is used
        :type h : float
        :return coordinates : [lon,lat,h] (3D np.array)
        """

        if h is not None:
            return self.model.direct_loc_h(row,col,h)
        elif self.dtm is not None:
            return self.model.direct_loc_dtm(row,col, self.dtm)
        else:
            return self.model.direct_loc_h(row, col, 0.0)

    def inverse(self,lon,lat,h):
        """
        inverse localization
        :param lat :  latitude
        :param lon : longitude
        :param h : altitude
        :return coordinates : [row,col,valid] (2D np.array), valid == 1 if coordinates is valid
        :rtype numpy.array
        """

        if self.use_rpc == False and not hasattr(self.model, 'pred_ofset_scale_lon'):
            self.model.estimate_inverse_loc_predictor()
        return self.model.inverse_loc(lon,lat,h)


def coloc(model1, model2, row, col, alt):
    """
    Colocalization : direct localization with model1, then inverse localization with model2

    :param model1: geometric model 1
    :type model1: shareloc.grid or  shareloc.rpc
    :param model2: geometric model 2
    :type model2: shareloc.grid or  shareloc.rpc
    :param row: sensor row
    :type row: int or 1D numpy array
    :param col: sensor col
    :type col: int or 1D numpy array
    :param alt: altitude
    :type alt: int or 1D numpy array
    :return: Corresponding sensor position [row, col, True] in the geometric model 2
    :rtype : Tuple(1D np.array row position, 1D np.array col position, 1D np.array True)
    """
    geometric_model1 = Localization(model1)
    geometric_model2 = Localization(model2)

    if not isinstance(row, (list, np.ndarray)):
        row = np.array([row])
        col = np.array([col])
        alt = np.array([alt])

    # Estimate longitudes, latitudes, altitudes using direct localization with model1
    ground_coord = geometric_model1.direct(row, col, alt)

    if len(ground_coord.shape) == 1:
        ground_coord = np.expand_dims(ground_coord, 0)

    # Estimate sensor position (row, col, altitude) using inverse localization with model2
    sensor_coord = np.zeros((row.shape[0], 3), dtype=np.float32)
    sensor_coord[:, 0], sensor_coord[:, 1], sensor_coord[:, 2] = geometric_model2.inverse(ground_coord[:, 0],
                                                                                          ground_coord[:, 1],
                                                                                          ground_coord[:, 2])

    return sensor_coord[:, 0], sensor_coord[:, 1], sensor_coord[:, 2]
