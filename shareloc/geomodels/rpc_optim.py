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
This module contains the optimized (with cpp bindings) RPC class corresponding to the RPC models.
RPC models covered are : DIMAP V1, DIMAP V2, DIMAP V3, ossim (geom file), geotiff.
"""


# Third party imports
from numba import config

import rpc_c

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.geomodels.rpc_readers import rpc_reader

# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"


@GeoModel.register("RpcOptim")
class RpcOptim(rpc_c.RPC, GeoModelTemplate):
    """
    RPC optimized with cpp bindings class including direct and inverse localization instance methods
    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self, rpc_params):
        GeoModelTemplate.__init__(self)

        norm_coeffs = [
            rpc_params["offset_x"],
            rpc_params["scale_x"],  # longitude
            rpc_params["offset_y"],
            rpc_params["scale_y"],  # latitude
            rpc_params["offset_alt"],
            rpc_params["scale_alt"],
            rpc_params["offset_col"],
            rpc_params["scale_col"],
            rpc_params["offset_row"],
            rpc_params["scale_row"],
        ]

        rpc_c.RPC.__init__(
            self,
            rpc_params["num_col"],
            rpc_params["den_col"],
            rpc_params["num_row"],
            rpc_params["den_row"],
            norm_coeffs,
        )

    @classmethod
    def load(cls, geomodel_path):
        """
        Load from any RPC (auto identify driver)
        from filename (dimap, ossim kwl, geotiff)

        TODO: topleftconvention always to True, set a standard and remove the option

        topleftconvention boolean: [0,0] position
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        # Set topleftconvention (keeping historic option): to clean
        cls.geomodel_path = geomodel_path
        return cls(rpc_reader(geomodel_path, topleftconvention=True))

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

    def inverse_loc(self, lon, lat, alt):
        """
        Inverse localization using c++ bindings

        :param lon: longitude position
        :type lon: float or 1D numpy.ndarray dtype=float64
        :param lat: latitude position
        :type lat: float or 1D numpy.ndarray dtype=float64
        :param alt: altitude
        :type alt: float
        :return: sensor position (row, col, alt)
        :rtype: tuple(1D np.array row position, 1D np.array col position, 1D np.array alt)
        """
        (row, col, alt) = super().inverse_loc(lon, lat, alt)
        return row, col, alt
