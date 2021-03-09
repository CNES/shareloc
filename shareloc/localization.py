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

"""
Localization class for localization functions.
"""

import numbers
import numpy as np


class Localization:
    """base class for localization function.
    Underlying model can be both multi layer localization grids or RPCs models
    """

    def __init__(self, model, elevation=None, image=None):
        """
        constructor
        :param model : geometric model
        :type model  : shareloc.grid or  shareloc.rpc
        :param elevation  : dtm or default elevation over ellipsoid if None elevation is set to 0
        :type elevation  : shareloc.dtm or float or np.ndarray
        :param image  : image class to handle geotransform
        :type image  : shareloc.image.image.Image
        """
        self.use_rpc = model.type == "rpc"
        self.model = model
        self.default_elevation = 0.0
        self.dtm = None
        if isinstance(elevation, (numbers.Number, list, np.ndarray)):
            self.default_elevation = elevation
        else:
            self.dtm = elevation
        self.image = image

    def direct(self, row, col, h=None, using_geotransform=False):
        """
        direct localization
        :param row :  sensor row
        :type row : float
        :param col : sensor col
        :type col : float
        :param h: altitude, if none DTM is used
        :type h : float
        :param using_geotransform: using_geotransform
        :type using_geotransform : boolean
        :return coordinates : [lon,lat,h] (3D np.array)
        """
        if using_geotransform and self.image is not None:
            row, col = self.image.transform_index_to_physical_point(row, col)
        if h is not None:
            return self.model.direct_loc_h(row, col, h)
        if self.dtm is not None:
            return self.model.direct_loc_dtm(row, col, self.dtm)
        return self.model.direct_loc_h(row, col, self.default_elevation)

    def inverse(self, lon, lat, h, using_geotransform=False):
        """
        inverse localization
        :param lat :  latitude
        :param lon : longitude
        :param h : altitude
        :param using_geotransform: using_geotransform
        :type using_geotransform : boolean
        :return coordinates : [row,col,h] (2D np.array)
        :rtype numpy.array
        """

        if not self.use_rpc and not hasattr(self.model, "pred_ofset_scale_lon"):
            self.model.estimate_inverse_loc_predictor()
        row, col, __ = self.model.inverse_loc(lon, lat, h)

        if using_geotransform and self.image is not None:
            row, col = self.image.transform_physical_point_to_index(row, col)
        return row, col, h


def coloc(model1, model2, row, col, elevation=None, image1=None, image2=None, using_geotransform=False):
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
    :param elevation: elevation
    :type elevation: shareloc.dtm or float or 1D numpy array
    :param image1  : image class to handle geotransform
    :type image1  : shareloc.image.image.Image
    :param image2  : image class to handle geotransform
    :type image2  : shareloc.image.image.Image
    :param using_geotransform: using_geotransform
    :type using_geotransform : boolean
    :return: Corresponding sensor position [row, col, True] in the geometric model 2
    :rtype : Tuple(1D np.array row position, 1D np.array col position, 1D np.array True)
    """
    geometric_model1 = Localization(model1, elevation, image=image1)
    geometric_model2 = Localization(model2, elevation, image=image2)

    if not isinstance(row, (list, np.ndarray)):
        row = np.array([row])
        col = np.array([col])

    ground_coord = geometric_model1.direct(row, col, using_geotransform=using_geotransform)

    if len(ground_coord.shape) == 1:
        ground_coord = np.expand_dims(ground_coord, 0)

    # Estimate sensor position (row, col, altitude) using inverse localization with model2
    sensor_coord = np.zeros((row.shape[0], 3), dtype=np.float64)
    sensor_coord[:, 0], sensor_coord[:, 1], sensor_coord[:, 2] = geometric_model2.inverse(
        ground_coord[:, 0], ground_coord[:, 1], ground_coord[:, 2], using_geotransform
    )

    return sensor_coord[:, 0], sensor_coord[:, 1], sensor_coord[:, 2]
