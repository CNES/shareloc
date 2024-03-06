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
This module contains the RPC class corresponding to the RPC models.
RPC models covered are : DIMAP V1, DIMAP V2, DIMAP V3, ossim (geom file), geotiff.
"""

# Standard imports
import logging
import os
from ast import literal_eval
from typing import Union

# Third party imports
import numpy as np
from affine import Affine
from numba import config, njit, prange

import bindings_cpp
from shareloc.geofunctions.dtm_intersection import DTMIntersection

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.geomodels.rpc_readers import rpc_reader
from shareloc.proj_utils import coordinates_conversion, transform_index_to_physical_point

# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"


@GeoModel.register("RPC")
class RPC(GeoModelTemplate):
    """
    RPC class including direct and inverse localization instance methods
    """

    # gitlab issue #61
    # pylint: disable=too-many-instance-attributes
    def __init__(self, rpc_params):
        super().__init__()

        self.offset_alt = None
        self.scale_alt = None
        self.offset_col = None
        self.scale_col = None
        self.offset_row = None
        self.scale_row = None
        self.offset_x = None
        self.scale_x = None
        self.offset_y = None
        self.scale_y = None

        self.datum = None
        for key, value in rpc_params.items():
            setattr(self, key, value)

        self.type = "RPC"
        if self.epsg is None:
            self.epsg = 4326
        if self.datum is None:
            self.datum = "ellipsoid"

        self.lim_extrapol = 1.0001

        # Each monome: c[0]*X**c[1]*Y**c[2]*Z**c[3]
        monomes_order = [
            [1, 0, 0, 0],
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [1, 0, 0, 1],
            [1, 1, 1, 0],
            [1, 1, 0, 1],
            [1, 0, 1, 1],
            [1, 2, 0, 0],
            [1, 0, 2, 0],
            [1, 0, 0, 2],
            [1, 1, 1, 1],
            [1, 3, 0, 0],
            [1, 1, 2, 0],
            [1, 1, 0, 2],
            [1, 2, 1, 0],
            [1, 0, 3, 0],
            [1, 0, 1, 2],
            [1, 2, 0, 1],
            [1, 0, 2, 1],
            [1, 0, 0, 3],
        ]

        self.monomes = np.array(monomes_order)

        # monomial coefficients of 1st variable derivative
        self.monomes_deriv_1 = np.array(
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [1, 0, 1, 0],
                [1, 0, 0, 1],
                [0, 0, 1, 1],
                [2, 1, 0, 0],
                [0, 0, 2, 0],
                [0, 0, 0, 2],
                [1, 0, 1, 1],
                [3, 2, 0, 0],
                [1, 0, 2, 0],
                [1, 0, 0, 2],
                [2, 1, 1, 0],
                [0, 0, 3, 0],
                [0, 0, 1, 2],
                [2, 1, 0, 1],
                [0, 0, 2, 1],
                [0, 0, 0, 3],
            ]
        )

        # monomial coefficients of 2nd variable derivative
        self.monomes_deriv_2 = np.array(
            [
                [0, 0, 0, 0],
                [0, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 0, 1],
                [1, 1, 0, 0],
                [0, 1, 0, 1],
                [1, 0, 0, 1],
                [0, 2, 0, 0],
                [2, 0, 1, 0],
                [0, 0, 0, 2],
                [1, 1, 0, 1],
                [0, 3, 0, 0],
                [2, 1, 1, 0],
                [0, 1, 0, 2],
                [1, 2, 0, 0],
                [3, 0, 2, 0],
                [1, 0, 0, 2],
                [0, 2, 0, 1],
                [2, 0, 1, 1],
                [0, 0, 0, 3],
            ]
        )

        self.inverse_coefficient = False
        self.direct_coefficient = False

        # pylint: disable=access-member-before-definition
        if self.num_col:
            self.inverse_coefficient = True
            self.num_col = np.array(self.num_col)
            self.den_col = np.array(self.den_col)
            self.num_row = np.array(self.num_row)
            self.den_row = np.array(self.den_row)

        # pylint: disable=access-member-before-definition
        if self.num_x:
            self.direct_coefficient = True
            self.num_x = np.array(self.num_x)
            self.den_x = np.array(self.den_x)
            self.num_y = np.array(self.num_y)
            self.den_y = np.array(self.den_y)

        self.alt_minmax = [self.offset_alt - self.scale_alt, self.offset_alt + self.scale_alt]
        self.col0 = self.offset_col - self.scale_col
        self.colmax = self.offset_col + self.scale_col
        self.row0 = self.offset_row - self.scale_row
        self.rowmax = self.offset_row + self.scale_row

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
        :param using_direct_coef: equals True if you want to use direct coefficients
        :type using_direct_coef: boolean
        :return: ground position (lon,lat,h)
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """
        if not isinstance(col, (list, np.ndarray)):
            col = np.array([col])
            row = np.array([row])

        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != col.shape[0]:
            alt = np.full(col.shape[0], fill_value=alt[0])

        points = np.zeros((col.size, 3))
        filter_nan, points[:, 0], points[:, 1] = self.filter_coordinates(row, col, fill_nan)
        row = row[filter_nan]
        col = col[filter_nan]

        # Direct localization using inverse RPC
        if not using_direct_coef and self.inverse_coefficient:
            logging.debug("direct localisation from inverse iterative")
            (points[filter_nan, 0], points[filter_nan, 1], points[filter_nan, 2]) = self.direct_loc_inverse_iterative(
                row, col, alt[filter_nan], 10, fill_nan
            )
        # Direct localization using direct RPC
        elif using_direct_coef and self.direct_coefficient:
            # ground position
            col_norm = (col - self.offset_col) / self.scale_col
            row_norm = (row - self.offset_row) / self.scale_row
            alt_norm = (alt[filter_nan] - self.offset_alt) / self.scale_alt

            if np.sum(abs(col_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! column extrapolation in direct localization ")
            if np.sum(abs(row_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! row extrapolation in direct localization ")
            if np.sum(abs(alt_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! alt extrapolation in direct localization ")

            points[filter_nan, 1], points[filter_nan, 0] = compute_rational_function_polynomial(
                col_norm,
                row_norm,
                alt_norm,
                self.num_x,
                self.den_x,
                self.num_y,
                self.den_y,
                self.scale_x,
                self.offset_x,
                self.scale_y,
                self.offset_y,
            )

        else:
            raise ValueError("Direct_loc_h: using_direct_coef doesn't match with available coefficients")

        points[:, 2] = alt
        return points

    def direct_loc_grid_h(self, row0, col0, steprow, stepcol, nbrow, nbcol, alt):
        """
        calculates a direct loc grid (lat, lon) from the direct RPCs at constant altitude
        TODO: not tested.

        :param row0:  grid origin (row)
        :type row0: int
        :param col0:  grid origin (col)
        :type col0: int
        :param steprow:  grid step (row)
        :type steprow: int
        :param stepcol:  grid step (col)
        :type stepcol: int
        :param nbrow:  grid nb row
        :type nbrow: int
        :param nbcol:  grid nb col
        :type nbcol: int
        :param alt: altitude of the grid
        :type alt: float
        :return: direct localization grid longitude and latitude
        :rtype: Tuple(numpy.array, numpy.array)
        """
        gri_lon = np.zeros((nbrow, nbcol))
        gri_lat = np.zeros((nbrow, nbcol))
        for column in range(int(nbcol)):
            col = col0 + stepcol * column
            for line in range(int(nbrow)):
                row = row0 + steprow * line
                (gri_lon[line, column], gri_lat[line, column], __) = self.direct_loc_h(row, col, alt)
        return (gri_lon, gri_lat)

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
        if not isinstance(col, (list, np.ndarray)):
            row = np.array([row])
            col = np.array([col])

        diff_alti_min, diff_alti_max = self.get_dtm_alt_offset(dtm.get_footprint_corners(), dtm)

        # print("min {} max {}".format(dtm.Zmin,dtm.Zmax))
        (min_dtm, max_dtm) = (dtm.get_alt_min() - 1.0 + diff_alti_min, dtm.get_alt_max() + 1.0 + diff_alti_max)
        if min_dtm < self.offset_alt - self.scale_alt:
            logging.debug("minimum dtm value is outside RPC validity domain, extrapolation will be done")
        if max_dtm > self.offset_alt + self.scale_alt:
            logging.debug("maximum dtm value is outside RPC validity domain, extrapolation will be done")

        los = self.los_extrema(row, col, min_dtm, max_dtm, epsg=dtm.get_epsg())

        # los -> (nb_point,nb_alt,3)
        los = los.T
        los = np.expand_dims(los, axis=0)
        los = np.moveaxis(los, 1, -1)
        los = los.reshape((len(col), 2, 3))
        direct_dtm = dtm.intersection_n_los_dtm(los)

        return direct_dtm

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
        if self.inverse_coefficient:
            if not isinstance(lon, (list, np.ndarray)):
                lon = np.array([lon])
                lat = np.array([lat])
            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != lon.shape[0]:
                alt = np.full(lon.shape[0], fill_value=alt[0])

            lon_norm = (lon - self.offset_x) / self.scale_x
            lat_norm = (lat - self.offset_y) / self.scale_y
            alt_norm = (alt - self.offset_alt) / self.scale_alt

            if np.sum(abs(lon_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! longitude extrapolation in inverse localization ")
            if np.sum(abs(lat_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! row extrapolation in inverse localization ")
            if np.sum(abs(alt_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! alt extrapolation in inverse localization ")

            row_out, col_out = compute_rational_function_polynomial(
                lon_norm,
                lat_norm,
                alt_norm,
                self.num_col,
                self.den_col,
                self.num_row,
                self.den_row,
                self.scale_col,
                self.offset_col,
                self.scale_row,
                self.offset_row,
            )
        else:
            logging.warning("inverse localisation can't be performed, inverse coefficients have not been defined")
            (col_out, row_out) = (None, None)
        return row_out, col_out, alt

    def filter_coordinates(self, first_coord, second_coord, fill_nan=False, direction="direct"):
        """
        Filter nan input values

        :param first_coord:  first coordinate
        :type first_coord: 1D numpy.ndarray dtype=float64
        :param second_coord:  second coordinate
        :type second_coord: 1D numpy.ndarray dtype=float64
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan: boolean
        :param direction:  direct or inverse localisation
        :type direction: str in ('direct', 'inverse')
        :return: filtered coordinates
        :rtype: list of numpy.array (index of nan, first filtered, second filtered)
        """
        filter_nan = np.logical_not(np.logical_or(np.isnan(first_coord), np.isnan(second_coord)))

        if fill_nan:
            if direction == "direct":
                out_x_nan_value = self.offset_x
                out_y_nan_value = self.offset_y
            else:
                out_x_nan_value = self.offset_col
                out_y_nan_value = self.offset_row
        else:
            out_x_nan_value = np.nan
            out_y_nan_value = np.nan

        x_out = np.full(len(second_coord), out_x_nan_value)
        y_out = np.full(len(second_coord), out_y_nan_value)

        return filter_nan, x_out, y_out

    def compute_loc_inverse_derivates(self, lon, lat, alt):
        """
        Inverse loc partial derivatives analytical compute

        :param lon: longitude coordinate
        :param lat: latitude coordinate
        :param alt: altitude coordinate
        :return: partials derivatives of inverse localization
        :rtype: Tuple(dcol_dlon np.array, dcol_dlat np.array, drow_dlon np.array, drow_dlat np.array)
        """
        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != lon.shape[0]:
            alt = np.full(lon.shape[0], fill_value=alt[0])

        lon_norm = (lon - self.offset_x) / self.scale_x
        lat_norm = (lat - self.offset_y) / self.scale_y
        alt_norm = (alt - self.offset_alt) / self.scale_alt

        dcol_dlon, dcol_dlat, drow_dlon, drow_dlat = compute_loc_inverse_derivates_numba(
            lon_norm,
            lat_norm,
            alt_norm,
            self.num_col,
            self.den_col,
            self.num_row,
            self.den_row,
            self.scale_col,
            self.scale_x,
            self.scale_row,
            self.scale_y,
        )

        return (dcol_dlon, dcol_dlat, drow_dlon, drow_dlat)

    def direct_loc_inverse_iterative(self, row, col, alt, nb_iter_max=10, fill_nan=False):
        """
        Iterative direct localization using inverse RPC

        :param row:  line sensor position
        :type row: float or 1D numpy.ndarray dtype=float64
        :param col:  column sensor position
        :type col: float or 1D numpy.ndarray dtype=float64
        :param alt:  altitude
        :type alt: float
        :param nb_iter_max: max number of iteration
        :type alt: int
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan: boolean
        :return: ground position (lon,lat,h)
        :rtype: list of numpy.array
        """

        if self.inverse_coefficient:
            if not isinstance(row, (list, np.ndarray)):
                col = np.array([col])
                row = np.array([row])

            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != col.shape[0]:
                alt = np.full(col.shape[0], fill_value=alt[0])

            filter_nan, long_out, lat_out = self.filter_coordinates(row, col, fill_nan)

            # if all coord contains Nan then return
            if not np.any(filter_nan):
                return long_out, lat_out, alt[filter_nan]

            row = row[filter_nan]
            col = col[filter_nan]
            alt_filtered = alt[filter_nan]

            # inverse localization starting from the center of the scene
            lon = np.array([self.offset_x])
            lat = np.array([self.offset_y])
            (row_start, col_start, __) = self.inverse_loc(lon, lat, alt_filtered)

            # desired precision in pixels
            eps = 1e-6

            iteration = 0
            # computing the residue between the sensor positions and those estimated by the inverse localization
            delta_col = col - col_start
            delta_row = row - row_start

            # ground coordinates (latitude and longitude) of each point
            lon = np.repeat(lon, delta_col.size)
            lat = np.repeat(lat, delta_col.size)

            # while the required precision is not achieved
            while (np.max(abs(delta_col)) > eps or np.max(abs(delta_row)) > eps) and iteration < nb_iter_max:
                # list of points that require another iteration
                iter_ = np.where((abs(delta_col) > eps) | (abs(delta_row) > eps))[0]

                # partial derivatives
                (dcol_dlon, dcol_dlat, drow_dlon, drow_dlat) = self.compute_loc_inverse_derivates(
                    lon[iter_], lat[iter_], alt_filtered[iter_]
                )

                det = dcol_dlon * drow_dlat - drow_dlon * dcol_dlat

                delta_lon = (drow_dlat * delta_col[iter_] - dcol_dlat * delta_row[iter_]) / det
                delta_lat = (-drow_dlon * delta_col[iter_] + dcol_dlon * delta_row[iter_]) / det

                # update ground coordinates
                lon[iter_] += delta_lon
                lat[iter_] += delta_lat

                # inverse localization
                (row_estim, col_estim, __) = self.inverse_loc(lon[iter_], lat[iter_], alt_filtered[iter_])

                # updating the residue between the sensor positions and those estimated by the inverse localization
                delta_col[iter_] = col[iter_] - col_estim
                delta_row[iter_] = row[iter_] - row_estim

                iteration += 1

            long_out[filter_nan] = lon
            lat_out[filter_nan] = lat

        else:
            logging.warning("inverse localisation can't be performed, inverse coefficients have not been defined")
            (long_out, lat_out) = (None, None)

        return long_out, lat_out, alt

    def get_alt_min_max(self):
        """
        returns altitudes min and max layers

        :return: alt_min,lat_max
        :rtype: list
        """
        return [self.offset_alt - self.scale_alt / 2.0, self.offset_alt + self.scale_alt / 2.0]

    def get_dtm_alt_offset(self, corners: np.ndarray, dtm: Union[DTMIntersection, bindings_cpp.DTMIntersection]):
        """
        returns min/max altitude offset between dtm coordinates system and RPC one

        :param corners: corners of the DTM's footprint
        :type corners: np.ndarray (4x2)
        :param dtm: DTM to get alt offset from
        :type dtm: DTMIntersection or bindings_cpp.DTMIntersection
        :return: min/max altimetric difference between RPC's epsg minus dtm alti expressed in dtm epsg
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

    def los_extrema(self, row, col, alt_min=None, alt_max=None, fill_nan=False, epsg=None):
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
        extrapolate = False
        if alt_min is None or alt_max is None:
            [los_alt_min, los_alt_max] = self.get_alt_min_max()
        elif alt_min >= self.alt_minmax[0] and alt_max <= self.alt_minmax[1]:
            los_alt_min = alt_min
            los_alt_max = alt_max
        else:
            extrapolate = True
            [los_alt_min, los_alt_max] = self.get_alt_min_max()

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
        los_edges = self.direct_loc_h(row_array, col_array, alt_array, fill_nan)

        if extrapolate:
            diff = los_edges[0::2, :] - los_edges[1::2, :]
            delta_alt = diff[:, 2]
            coeff_alt_max = (alt_max - los_edges[1::2, 2]) / delta_alt
            coeff_alt_max = np.tile(coeff_alt_max[:, np.newaxis], (1, 3))
            coeff_alt_min = (alt_min - los_edges[1::2, 2]) / delta_alt
            coeff_alt_min = np.tile(coeff_alt_min[:, np.newaxis], (1, 3))
            los_edges[0::2, :] = los_edges[1::2, :] + diff * coeff_alt_max
            los_edges[1::2, :] = los_edges[1::2, :] + diff * coeff_alt_min

        if epsg is not None and epsg != self.epsg:
            los_edges = coordinates_conversion(los_edges, self.epsg, epsg)

        return los_edges


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def polynomial_equation(xnorm, ynorm, znorm, coeff):
    """
    Compute polynomial equation

    :param xnorm: Normalized longitude (for inverse) or column (for direct) position
    :type xnorm: float 64
    :param ynorm: Normalized latitude (for inverse) or line (for direct) position
    :type ynorm: float 64
    :param znorm: Normalized altitude position
    :type znorm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational
    :rtype: float 64
    """
    rational = (
        coeff[0]
        + coeff[1] * xnorm
        + coeff[2] * ynorm
        + coeff[3] * znorm
        + coeff[4] * xnorm * ynorm
        + coeff[5] * xnorm * znorm
        + coeff[6] * ynorm * znorm
        + coeff[7] * xnorm**2
        + coeff[8] * ynorm**2
        + coeff[9] * znorm**2
        + coeff[10] * xnorm * ynorm * znorm
        + coeff[11] * xnorm**3
        + coeff[12] * xnorm * ynorm**2
        + coeff[13] * xnorm * znorm**2
        + coeff[14] * xnorm**2 * ynorm
        + coeff[15] * ynorm**3
        + coeff[16] * ynorm * znorm**2
        + coeff[17] * xnorm**2 * znorm
        + coeff[18] * ynorm**2 * znorm
        + coeff[19] * znorm**3
    )

    return rational


# pylint: disable=too-many-arguments
@njit(
    "Tuple((f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)",
    parallel=literal_eval(os.environ.get("SHARELOC_NUMBA_PARALLEL", "True")),
    cache=True,
    fastmath=True,
)
def compute_rational_function_polynomial(
    lon_col_norm,
    lat_row_norm,
    alt_norm,
    num_col,
    den_col,
    num_lin,
    den_lin,
    scale_col,
    offset_col,
    scale_lin,
    offset_lin,
):
    """
    Compute rational function polynomial using numba to reduce calculation time on multiple points.
    useful to compute direct and inverse localization using direct or inverse RPC.

    :param lon_col_norm: Normalized longitude (for inverse) or column (for direct) position
    :type lon_col_norm: 1D np.array dtype np.float 64
    :param lat_row_norm: Normalized latitude (for inverse) or line (for direct) position
    :type lat_row_norm: 1D np.array dtype np.float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: 1D np.array dtype np.float 64
    :param num_col: Column numerator coefficients
    :type num_col: 1D np.array dtype np.float 64
    :param den_col: Column denominator coefficients
    :type den_col: 1D np.array dtype np.float 64
    :param num_lin: Line numerator coefficients
    :type num_lin: 1D np.array dtype np.float 64
    :param den_lin: Line denominator coefficients
    :type den_lin: 1D np.array dtype np.float 64
    :param scale_col: Column scale
    :type scale_col: float 64
    :param offset_col: Column offset
    :type offset_col: float 64
    :param scale_lin: Line scale
    :type scale_lin: float 64
    :param offset_lin: Line offset
    :type offset_lin: float 64
    :return: for inverse localization : sensor position (row, col). for direct localization : ground position (lon, lat)
    :rtype: Tuple(np.ndarray, np.ndarray)
    """
    assert lon_col_norm.shape == alt_norm.shape

    col_lat_out = np.zeros((lon_col_norm.shape[0]), dtype=np.float64)
    row_lon_out = np.zeros((lon_col_norm.shape[0]), dtype=np.float64)

    # pylint: disable=not-an-iterable
    for i in prange(lon_col_norm.shape[0]):
        poly_num_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_col)
        poly_den_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_col)
        poly_num_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_lin)
        poly_den_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_lin)
        col_lat_out[i] = poly_num_col / poly_den_col * scale_col + offset_col
        row_lon_out[i] = poly_num_lin / poly_den_lin * scale_lin + offset_lin

    return row_lon_out, col_lat_out


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, coeff):
    """
    Compute latitude derivative polynomial equation

    :param lon_norm: Normalized longitude position
    :type lon_norm: float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    derivate = (
        coeff[2]
        + coeff[4] * lon_norm
        + coeff[6] * alt_norm
        + 2 * coeff[8] * lat_norm
        + coeff[10] * lon_norm * alt_norm
        + 2 * coeff[12] * lon_norm * lat_norm
        + coeff[14] * lon_norm**2
        + 3 * coeff[15] * lat_norm**2
        + coeff[16] * alt_norm**2
        + 2 * coeff[18] * lat_norm * alt_norm
    )

    return derivate


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, coeff):
    """
    Compute longitude derivative polynomial equation

    :param lon_norm: Normalized longitude position
    :type lon_norm: float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    derivate = (
        coeff[1]
        + coeff[4] * lat_norm
        + coeff[5] * alt_norm
        + 2 * coeff[7] * lon_norm
        + coeff[10] * lat_norm * alt_norm
        + 3 * coeff[11] * lon_norm**2
        + coeff[12] * lat_norm**2
        + coeff[13] * alt_norm**2
        + 2 * coeff[14] * lat_norm * lon_norm
        + 2 * coeff[17] * lon_norm * alt_norm
    )

    return derivate


# pylint: disable=too-many-arguments
@njit(
    "Tuple((f8[:], f8[:], f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)",
    parallel=literal_eval(os.environ.get("SHARELOC_NUMBA_PARALLEL", "True")),
    cache=True,
    fastmath=True,
)
def compute_loc_inverse_derivates_numba(
    lon_norm, lat_norm, alt_norm, num_col, den_col, num_lin, den_lin, scale_col, scale_lon, scale_lin, scale_lat
):
    """
    Analytically compute the partials derivatives of inverse localization using numba to reduce calculation time on
    multiple points

    :param lon_norm: Normalized longitude position
    :type lon_norm: 1D np.array dtype np.float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: 1D np.array dtype np.float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: 1D np.array dtype np.float 64
    :param num_col: Column numerator coefficients
    :type num_col: 1D np.array dtype np.float 64
    :param den_col: Column denominator coefficients
    :type den_col: 1D np.array dtype np.float 64
    :param num_lin: Line numerator coefficients
    :type num_lin: 1D np.array dtype np.float 64
    :param den_lin: Line denominator coefficients
    :type den_lin: 1D np.array dtype np.float 64
    :param scale_col: Column scale
    :type scale_col: float 64
    :param scale_lon: Geodetic longitude scale
    :type scale_lon: float 64
    :param scale_lin: Line scale
    :type scale_lin: float 64
    :param scale_lat: Geodetic latitude scale
    :type scale_lat: float 64
    :return: partials derivatives of inverse localization
    :rtype: Tuples(dcol_dlon np.array, dcol_dlat np.array, drow_dlon np.array, drow_dlat np.array)
    """
    dcol_dlon = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    dcol_dlat = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    drow_dlon = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    drow_dlat = np.zeros((lon_norm.shape[0]), dtype=np.float64)

    # pylint: disable=not-an-iterable
    for i in prange(lon_norm.shape[0]):
        num_dcol = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        num_dcol_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        num_dcol_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        dcol_dlon[i] = scale_col / scale_lon * (num_dcol_dlon * den_dcol - den_dcol_dlon * num_dcol) / den_dcol**2
        dcol_dlat[i] = scale_col / scale_lat * (num_dcol_dlat * den_dcol - den_dcol_dlat * num_dcol) / den_dcol**2
        drow_dlon[i] = scale_lin / scale_lon * (num_drow_dlon * den_drow - den_drow_dlon * num_drow) / den_drow**2
        drow_dlat[i] = scale_lin / scale_lat * (num_drow_dlat * den_drow - den_drow_dlat * num_drow) / den_drow**2

    return dcol_dlon, dcol_dlat, drow_dlon, drow_dlat
