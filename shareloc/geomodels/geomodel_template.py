#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of shareloc
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
This module contains the GeoModel abstract class
"""

# Standard imports
from abc import abstractmethod


class GeoModelTemplate:
    """
    Class for general specification of a geometric model
    declined in rpc.py and grid.py and rpc_optim.py
    """

    @abstractmethod
    def __init__(self):
        """
        Return the geomodel object associated with the geomodel_type
        given in the configuration
        """
        # geomodel type. Set by the subclass
        self.type: str
        # geomodel epsg projection code
        self.epsg: int = None

    # Define GeoModelTemplate functions interface

    @abstractmethod
    def direct_loc_h(self, row, col, alt, fill_nan=False):
        """
        direct localization at constant altitude

        :param row:  line sensor position
        :type row: float or 1D numpy.ndarray dtype=float64
        :param col:  column sensor position
        :type col: float or 1D numpy.ndarray dtype=float64
        :param alt:  altitude
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan: boolean
        :return: ground position (lon,lat,h)
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """

    @abstractmethod
    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm

        :param row:  line sensor position
        :type row: float
        :param col:  column sensor position
        :type col: float
        :param dtm: dtm intersection model
        :type dtm: shareloc.geofunctions.dtm_intersection
        :return: ground position (lon,lat,h) in dtm coordinates system
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """

    @abstractmethod
    def inverse_loc(self, lon, lat, alt):
        """
        Inverse localization

        :param lon: longitude position
        :type lon: float or 1D numpy.ndarray dtype=float64
        :param lat: latitude position
        :type lat: float or 1D numpy.ndarray dtype=float64
        :param alt: altitude
        :type alt: float
        :return: sensor position (row, col, alt)
        :rtype: tuple(1D np.array row position, 1D np.array col position, 1D np.array alt)
        """

    @classmethod
    @abstractmethod
    def load(cls, geomodel_path: str):
        """
        load function with class specific args

        :param geomodel_path: filename of geomodel
        """
