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


import logging
from typing import Union

import numpy as np
from affine import Affine

# Third party imports
from numba import config

import bindings_cpp
from shareloc.geofunctions.dtm_intersection import DTMIntersection

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.geomodels.rpc_readers import rpc_reader
from shareloc.proj_utils import coordinates_conversion, transform_index_to_physical_point

# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"


@GeoModel.register("RPCoptim")
class RPCoptim(bindings_cpp.RPC, GeoModelTemplate):
    """
    RPC optimized with cpp bindings class including direct and inverse localization instance methods
    """

    # pylint: disable=too-many-instance-attributes,duplicate-code
    def __init__(self, rpc_params):
        GeoModelTemplate.__init__(self)

        self.type = "RPCoptim"
        if "espg" not in rpc_params:
            self.epsg = 4326
        else:
            self.epsg = rpc_params["epsg"]

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

        empty = [0 for i in range(20)]

        if rpc_params["num_col"] and rpc_params["num_x"]:
            if not all(i == 0 for i in rpc_params["num_col"]) and not all(
                i == 0 for i in rpc_params["num_x"]
            ):  # direct and inverse coef
                bindings_cpp.RPC.__init__(
                    self,
                    True,
                    True,
                    rpc_params["num_col"],
                    rpc_params["den_col"],
                    rpc_params["num_row"],
                    rpc_params["den_row"],
                    rpc_params["num_x"],
                    rpc_params["den_x"],
                    rpc_params["num_y"],
                    rpc_params["den_y"],
                    norm_coeffs,
                )
            else:
                raise ValueError("RPCoptim : RPC coefficients are all 0")

        elif rpc_params["num_col"]:
            if not all(i == 0 for i in rpc_params["num_col"]):  # only inverse coef
                bindings_cpp.RPC.__init__(
                    self,
                    True,
                    False,
                    rpc_params["num_col"],
                    rpc_params["den_col"],
                    rpc_params["num_row"],
                    rpc_params["den_row"],
                    empty,
                    empty,
                    empty,
                    empty,
                    norm_coeffs,
                )
            else:
                raise ValueError("RPCoptim : RPC coefficients are all 0")

        elif rpc_params["num_x"]:
            if not all(i == 0 for i in rpc_params["num_x"]):  # only direct coef
                bindings_cpp.RPC.__init__(
                    self,
                    False,
                    True,
                    empty,
                    empty,
                    empty,
                    empty,
                    rpc_params["num_x"],
                    rpc_params["den_x"],
                    rpc_params["num_y"],
                    rpc_params["den_y"],
                    norm_coeffs,
                )
            else:
                raise ValueError("RPCoptim : RPC coefficients are all 0")
        else:
            raise ValueError("RPCoptim : No RPC coefficients readable")

        self.col0 = self.get_offset_col() - self.get_scale_col()
        self.colmax = self.get_offset_col() + self.get_scale_col()
        self.row0 = self.get_offset_row() - self.get_scale_row()
        self.rowmax = self.get_offset_row() + self.get_scale_row()

        self.epsg = 4326  # rpc_params does not have an espg key
        # But PHRDIMAP_ files have an epsg code specfied

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

    def direct_loc_h(self, row, col, alt, fill_nan=False, using_direct_coef=False):
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
        if not isinstance(alt, (list, np.ndarray)) and isinstance(row, (list, np.ndarray)):
            alt = [alt]

        # TODO : use py::array in c++
        res_optim = super().direct_loc_h(row, col, alt, fill_nan, using_direct_coef)
        # same output as python
        if not isinstance(res_optim[0], (list, np.ndarray)):
            res_optim_0 = [res_optim[0]]
            res_optim_1 = [res_optim[1]]
            res_optim_2 = [res_optim[2]]
        else:
            res_optim_0 = res_optim[0]
            res_optim_1 = res_optim[1]
            res_optim_2 = res_optim[2]

        res_optim = np.array([res_optim_0, res_optim_1, res_optim_2]).T

        return res_optim

    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm only if dtm's epsg code is 4326

        :param row:  line sensor position
        :type row: list or np.array
        :param col:  column sensor position
        :type col: list or np.array
        :param dtm: dtm intersection c++ model
        :type dtm: shareloc.bindings.dtm_intersection.cpp
        :return: ground position (lon,lat,h) in dtm coordinates system
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """

        if dtm.get_epsg() == 4326:  # full c++
            res_optim = super().direct_loc_dtm(row, col, dtm)

            if not isinstance(res_optim[0], (list, np.ndarray)):
                res_optim_0 = [res_optim[0]]
                res_optim_1 = [res_optim[1]]
                res_optim_2 = [res_optim[2]]
            else:
                res_optim_0 = res_optim[0]
                res_optim_1 = res_optim[1]
                res_optim_2 = res_optim[2]

            res_optim = np.array([res_optim_0, res_optim_1, res_optim_2]).T

        else:  # Beginning in python and core in c++

            diff_alti_min, diff_alti_max = self.get_dtm_alt_offset(dtm.get_footprint_corners(), dtm)

            if not isinstance(col, (list, np.ndarray)):
                row = np.array([row])
                col = np.array([col])

            (min_dtm, max_dtm) = (dtm.get_alt_min() - 1.0 + diff_alti_min, dtm.get_alt_max() + 1.0 + diff_alti_max)
            if min_dtm < self.get_offset_alt() - self.get_scale_alt():
                logging.debug("minimum dtm value is outside RPC validity domain, extrapolation will be done")
            if max_dtm > self.get_offset_alt() + self.get_scale_alt():
                logging.debug("maximum dtm value is outside RPC validity domain, extrapolation will be done")

            los = self.los_extrema(row, col, min_dtm, max_dtm, epsg=dtm.get_epsg())

            # los -> (nb_point,nb_alt,3)
            los = los.T
            los = np.expand_dims(los, axis=0)
            los = np.moveaxis(los, 1, -1)
            los = los.reshape((len(col), 2, 3))
            res_optim = dtm.intersection_n_los_dtm(los)

        return res_optim

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

        return np.array(row), np.array(col), np.array(alt)

    def get_dtm_alt_offset(self, corners: np.ndarray, dtm: Union[DTMIntersection, bindings_cpp.DTMIntersection]):
        """
        returns min/max altitude offset between dtm coordinates system and RPC one

        :param corners: corners of the DTM's footprint
        :type corners: np.ndarray (4x2)
        :param dtm: DTM to get alt offset from
        :type dtm: bindings.DTMIntersection
        :return: min/max altimetric difference between RPC'sepsg minus dtm alti expressed in dtm epsg
        :rtype: list of float (1x2)
        """

        alti_moy = (dtm.get_alt_min() + dtm.get_alt_max()) / 2.0

        ground_corners = np.zeros([3, 4])
        ground_corners[2, :] = alti_moy
        transform = dtm.get_transform()
        if not isinstance(transform, Affine):
            transform = Affine.from_gdal(*transform)
        ground_corners[1::-1, :] = transform_index_to_physical_point(transform, corners[:, 0], corners[:, 1])
        converted_corners = coordinates_conversion(ground_corners.transpose(), dtm.get_epsg(), self.epsg)
        return [np.min(converted_corners[:, 2]) - alti_moy, np.max(converted_corners[:, 2]) - alti_moy]

    def los_extrema(self, row, col, alt_min=None, alt_max=None, fill_nan=False, epsg=4326):
        """
        compute los extrema

        :param row:  line sensor position
        :type row: float
        :param col:  column sensor position
        :type col: float
        :param alt_min: los alt min
        :type alt_min: float
        :param alt_max: los alt max
        :type alt_max: float
        :param epsg: epsg code of the dtm
        :type epsg: int
        :return: los extrema
        :rtype: numpy.array (2x3)
        """
        if epsg == 4326:  # full c++
            res_cpp = super().los_extrema(row, col, alt_min, alt_max, fill_nan)
            res_cpp = np.array(res_cpp).T
            return res_cpp

        extrapolate = False
        if alt_min is None or alt_max is None:
            [los_alt_min, los_alt_max] = super().get_alt_min_max()
        elif alt_min >= super().get_alt_minmax()[0] and alt_max <= super().get_alt_minmax()[1]:
            los_alt_min = alt_min
            los_alt_max = alt_max
        else:
            extrapolate = True
            [los_alt_min, los_alt_max] = super().get_alt_min_max()

        if isinstance(row, (np.ndarray)):
            los_nb = row.shape[0]
            row_array = np.full([los_nb * 2], fill_value=0.0)
            col_array = np.full([los_nb * 2], fill_value=0.0)
            alt_array = np.full([los_nb * 2], fill_value=0.0)
            row_array[0::2] = row
            row_array[1::2] = row
            col_array[0::2] = col
            col_array[1::2] = col
            alt_array[0::2] = los_alt_max
            alt_array[1::2] = los_alt_min
        else:
            los_nb = 1
            row_array = np.array([row, row])
            col_array = np.array([col, col])
            alt_array = np.array([los_alt_max, los_alt_min])

        row_dh, col_dh, alt_dh = super().direct_loc_h(row_array, col_array, alt_array, fill_nan, False)
        los_edges = np.array([row_dh, col_dh, alt_dh]).T
        if extrapolate:
            diff = los_edges[0::2, :] - los_edges[1::2, :]
            delta_alt = diff[:, 2]
            coeff_alt_max = (alt_max - los_edges[1::2, 2]) / delta_alt
            coeff_alt_max = np.tile(coeff_alt_max[:, np.newaxis], (1, 3))
            coeff_alt_min = (alt_min - los_edges[1::2, 2]) / delta_alt
            coeff_alt_min = np.tile(coeff_alt_min[:, np.newaxis], (1, 3))
            los_edges[0::2, :] = los_edges[1::2, :] + diff * coeff_alt_max
            los_edges[1::2, :] = los_edges[1::2, :] + diff * coeff_alt_min

        los_edges = coordinates_conversion(los_edges, self.epsg, epsg)  # self.epsg = 4326

        return los_edges
