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
Localisation functions from multi h direct grids.
"""

# Standard imports
import logging

# Third party imports
import numpy as np

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.image import Image
from shareloc.math_utils import interpol_bilin_grid, interpol_bilin_vectorized
from shareloc.proj_utils import coordinates_conversion


# gitlab issue #58
# pylint: disable=too-many-instance-attributes
@GeoModel.register("GRID")
class Grid(GeoModelTemplate):
    """
    multi H direct localization grid handling class.
    please refer to the main documentation grid format

    Derives from GeoModelTemplate

    :param row0: grid first pixel center along Y axis (row).
    :type row0: float
    :param col0: grid first pixel center along X axis (column).
    :type col0: float
    :param nbrow: grid size in row
    :type nbrow: int
    :param nbcol: grid size in col
    :type nbcol: int
    :param steprow: grid step in row
    :type steprow: float
    :param stepcol: grid step in col
    :type stepcol: float
    :param rowmax: last row in grid
    :type rowmax: float
    :param colmax: last col in grid
    :type colmax: float
    :param repter: ground coordinate system
    :type repter: str
    :param epsg: epsg code corresponding to shareloc.grid.Grid.repter
    :type epsg: int
    :param nbalt: number of altitude layers
    :type nbalt: int
    :param lon_data: longitude array
    :type lon_data: np.ndarray of size (nbalt,nbrow,nbcol) size
    :param lat_data: latitude array
    :type lat_data: np.ndarray of size (nbalt,nbrow,nbcol) size
    :param alts_down: altitudes in decreasing order
    :type alts_down: list
    :param type: geometric model type
    :type type: str
    """

    def __init__(self, geomodel_path: str):
        """
        Grid Constructor

        :param geomodel_path: grid filename (Geotiff)
        :type geomodel_path: string
        """
        # Instanciate GeoModelTemplate generic init with shared parameters
        super().__init__()
        self.type = "multi H grid"
        self.geomodel_path = geomodel_path
        # GeoModel Grid parameters definition (see documentation)
        self.row0 = None
        self.col0 = None
        self.nbrow = None
        self.nbcol = None
        self.steprow = None
        self.stepcol = None
        self.repter = None
        self.nbalt = None
        self.lon_data = None  # longitude array
        self.lat_data = None  # latitude array
        self.alts_down = None
        self.rowmax = None
        self.colmax = None

        # inverse loc predictor attributes
        self.pred_col_min = None
        self.pred_row_min = None
        self.pred_col_max = None
        self.pred_row_max = None
        self.pred_ofset_scale_lon = None
        self.pred_ofset_scale_lat = None
        self.pred_ofset_scale_row = None
        self.pred_ofset_scale_col = None

        self.read()

    @classmethod
    def load(cls, geomodel_path):
        """
        Load grid and fill Class attributes.

        The grid is read as an shareloc.image. Image and class attributes are filled.
        Shareloc geotiff grids are stored by increasing altitude H0 ... Hx
        2 data cubes are defined:
        - lon_data : [alt,row,col]
        - lat_data : [alt,row,col]
        """
        return cls(geomodel_path)

    def read(self):
        """
        Load grid and fill Class attributes.

        The grid is read as an shareloc.image. Image and class attributes are filled.
        Shareloc geotiff grids are stored by increasing altitude H0 ... Hx
        2 data cubes are defined:
        - lon_data : [alt,row,col]
        - lat_data : [alt,row,col]
        """
        grid_image = Image(self.geomodel_path, read_data=True)
        if grid_image.dataset.driver != "GTiff":
            raise TypeError(
                "Only Geotiff grids are accepted. Please refer to the documentation for grid supported format."
            )
        metadata = grid_image.dataset.tags()
        self.nbalt = int(grid_image.dataset.count / 2)
        self.nbrow = grid_image.nb_rows
        self.nbcol = grid_image.nb_columns
        # data in global or specific metadata namespace
        for tags in grid_image.dataset.tag_namespaces():
            metadata.update(grid_image.dataset.tags(ns=tags))
        indexes = self.parse_metadata_alti(metadata)
        lon_indexes = indexes * 2
        self.lon_data = grid_image.data[lon_indexes, :, :]
        self.lat_data = grid_image.data[lon_indexes + 1, :, :]
        self.stepcol = grid_image.pixel_size_col
        self.steprow = grid_image.pixel_size_row
        self.col0 = grid_image.origin_col
        self.row0 = grid_image.origin_row
        self.rowmax = self.row0 + self.steprow * (self.nbrow - 1)
        self.colmax = self.col0 + self.stepcol * (self.nbcol - 1)
        for key in metadata:
            if key.endswith("REF"):
                proj_ref = metadata[key]
                self.repter = proj_ref
                if proj_ref.startswith("EPSG:"):
                    self.epsg = int(proj_ref.split(":")[1])
                else:
                    logging.debug("use default epsg : 4326 for the grid crs.")
                    self.epsg = 4326

    def parse_metadata_alti(self, metadata):
        """
        parse metadata to sort altitude in decreasing order

        :param metadata: Geotiff metadata
        :type metadata: dict
        """
        alt_array = np.zeros([self.nbalt])
        for key in metadata.keys():
            if "ALTITUDE" in key:
                band_index = int(int(key.split("B")[1]) / 2)
                alt_array[band_index] = float(metadata[key])
        indexes = np.argsort(alt_array)[::-1]
        self.alts_down = alt_array[indexes]
        return indexes

    def get_alt_min_max(self):
        """
        returns altitudes min and max layers

        :return: alt_min,lat_max
        :rtype: list
        """
        return [self.alts_down[-1], self.alts_down[0]]

    def direct_loc_h(self, row, col, alt, fill_nan=False):
        """
        direct localization at constant altitude

        :param row: line sensor position
        :type row: float or 1D numpy.ndarray dtype=float64
        :param col: column sensor position
        :type col: float or 1D numpy.ndarray dtype=float64
        :param alt: altitude
        :type alt: float
        :param fill_nan: not used, preserved for API symmetry
        :type fill_nan: boolean
        :return: ground position (lon,lat,h)
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """

        # Vectorization doesn't handle yet altitude as np.ndarray (3D interpolation work).
        # TODO: clean interfaces to remove this part
        if isinstance(alt, (list, np.ndarray)):
            logging.debug("grid doesn't handle alt as array, first value is used")
            alt = alt[0]
        if fill_nan:
            logging.warning("fill nan strategy not available for grids")
        (grid_index_up, grid_index_down) = self.return_grid_index(alt)
        alt_down = self.alts_down[grid_index_down]
        alt_up = self.alts_down[grid_index_up]
        alti_coef = (alt - alt_down) / (alt_up - alt_down)
        mats = [
            self.lon_data[grid_index_up : grid_index_down + 1, :, :],
            self.lat_data[grid_index_up : grid_index_down + 1, :, :],
        ]

        # float are converted to np.ndarray for vectorization
        # TODO: refactoring to remove this part.
        if not isinstance(col, (list, np.ndarray)):
            col = np.array([col])
            row = np.array([row])

        filter_nan = np.logical_not(np.logical_or(np.isnan(col), np.isnan(row)))
        position = np.nan * np.zeros((col.size, 3))
        position[:, 2] = alt
        if np.any(filter_nan):
            pos_row = (row[filter_nan] - self.row0) / self.steprow
            pos_col = (col[filter_nan] - self.col0) / self.stepcol
            # pylint disable for code clarity interpol_bilin_vectorized returns one list of 2 elements in this case
            # pylint: disable=unbalanced-tuple-unpacking
            [vlon, vlat] = interpol_bilin_vectorized(mats, self.nbrow, self.nbcol, pos_row, pos_col)
            position[filter_nan, 0] = alti_coef * vlon[0, :] + (1 - alti_coef) * vlon[1, :]
            position[filter_nan, 1] = alti_coef * vlat[0, :] + (1 - alti_coef) * vlat[1, :]
        return position

    def compute_los(self, row, col, epsg):
        """
        Compute Line of Sight

        :param row: line sensor position
        :type row: float
        :param col: column sensor position
        :type col: float
        :param epsg: epsg code
        :type epsg: int
        :return: los
        :rtype: numpy.array
        """
        los = np.zeros((3, self.nbalt))
        loslonlat = self.interpolate_grid_in_plani(row, col)
        los[0, :] = loslonlat[0]
        los[1, :] = loslonlat[1]
        los[2, :] = self.alts_down
        los = los.T
        if epsg != self.epsg:
            los = coordinates_conversion(los, self.epsg, epsg)
        return los

    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm

        TODO explain algorithm
        TODO optimize code (for loop, ...)

        :param row: line sensor position
        :type row: float
        :param col: column sensor position
        :type col: float
        :param dtm: dtm model
        :type dtm: shareloc.dtm
        :return: ground position (lon,lat,h) in dtm coordinates system.
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """
        if not isinstance(row, (list, np.ndarray)):
            row = np.array([row])
            col = np.array([col])

        filter_nan = np.logical_not(np.logical_or(np.isnan(col), np.isnan(row)))
        points_dtm = np.nan * np.zeros((col.size, 3))
        if np.any(filter_nan):
            row_filtered = row[filter_nan]
            col_filtered = col[filter_nan]
            points_nb = row_filtered.size
            all_los = np.empty((points_nb, self.nbalt, 3))
            for point_index in np.arange(points_nb):
                row_i = row_filtered[point_index]
                col_i = col_filtered[point_index]
                los = self.compute_los(row_i, col_i, dtm.get_epsg())

                all_los[point_index, :, :] = los

            points_dtm[filter_nan, :] = dtm.intersection_n_los_dtm(all_los)
        return points_dtm

    def los_extrema(self, row, col, alt_min, alt_max):
        """
        compute los extrema

        :param row: line sensor position
        :type row: float
        :param col: column sensor position
        :type col: float
        :param alt_min: los alt min
        :type alt_min: float
        :param alt_max: los alt max
        :type alt_max: float
        :return: los extrema
        :rtype: numpy.array (2x3)
        """
        los_edges = np.zeros([2, 3])
        los_edges[0, :] = self.direct_loc_h(row, col, alt_max)
        los_edges[1, :] = self.direct_loc_h(row, col, alt_min)
        return los_edges

    def interpolate_grid_in_plani(self, row, col):
        """
        interpolate positions on multi h grid

        :param row: line sensor position
        :type row: float
        :param col: column sensor position
        :type col: float
        :return: interpolated positions
        :rtype: list
        """
        pos_row = (row - self.row0) / self.steprow
        pos_col = (col - self.col0) / self.stepcol
        mats = [self.lon_data, self.lat_data]
        res = interpol_bilin_grid(mats, self.nbrow, self.nbcol, pos_row, pos_col)
        return res

    def interpolate_grid_in_altitude(self, nbrow, nbcol, nbalt=None):
        """
        interpolate equally spaced grid (in altitude)

        :param nbrow: grid nb row
        :type nbrow: int
        :param nbcol: grid nb col
        :type nbcol: int
        :param nbalt: grid nb alt, of None self.nbalt is used instead
        :type nbalt: int
        :return: equally spaced grid
        :rtype: numpy.array
        """
        if not nbalt:
            nbalt = self.nbalt
            list_alts = self.alts_down
        else:
            list_alts = np.linspace(self.alts_down[0], self.alts_down[-1], nbalt)

        lon_data = np.zeros((nbalt, nbrow, nbcol))
        lat_data = np.zeros((nbalt, nbrow, nbcol))

        # generates an interpolated direction cube of nrow/ncol directions
        # row_max = self.row0 + self.steprow * (self.nbrow-1)
        # col_max = self.col0 + self.stepcol * (self.nbcol-1)

        steprow = (self.rowmax - self.row0) / (nbrow - 1)
        stepcol = (self.colmax - self.col0) / (nbcol - 1)

        for index, alt in enumerate(list_alts):
            res = self.direct_loc_grid_h(self.row0, self.col0, steprow, stepcol, nbrow, nbcol, alt)
            lon_data[index] = res[0]
            lat_data[index] = res[1]
        return lon_data, lat_data

    def direct_loc_grid_dtm(self, row0, col0, steprow, stepcol, nbrow, nbcol, dtm):
        """
        direct localization  grid on dtm

        :param row0: grid origin (row)
        :type row0: int
        :param col0: grid origin (col)
        :type col0: int
        :param steprow: grid step (row)
        :type steprow: int
        :param stepcol: grid step (col)
        :type stepcol: int
        :param nbrow: grid nb row
        :type nbrow: int
        :param nbcol: grid nb col
        :type nbcol: int
        :param dtm: dtm model
        :type dtm: shareloc.dtm
        :return: direct localization grid
        :rtype: numpy.array
        """
        glddtm = np.zeros((3, nbrow, nbcol))
        los = np.zeros((3, self.nbalt))
        for i in range(nbrow):
            for j in range(nbcol):
                col = col0 + stepcol * j
                row = row0 + steprow * i
                los = self.compute_los(row, col, dtm.get_epsg())
                (__, position_cube, alti, los_index) = dtm.intersect_dtm_cube(los)
                if position_cube is not None:
                    (__, point_dtm) = dtm.intersection(los_index, position_cube, alti)
                else:
                    point_dtm = np.full(3, fill_value=np.nan)
                # conversion of all tab
                # if self.epsg != dtm.epsg:
                #    point_r = coordinates_conversion(point_r, dtm.epsg, self.epsg)
                glddtm[:, i, j] = point_dtm
        return glddtm

    def return_grid_index(self, alt):
        """
        return layer index enclosing a given altitude

        :param alt: altitude
        :type alt: float
        :return: grid index (up,down)
        :rtype: tuple
        """
        if alt > self.alts_down[0]:
            (high_index, low_index) = (0, 0)
        elif alt < self.alts_down[-1]:
            (high_index, low_index) = (self.nbalt - 1, self.nbalt - 1)
        else:
            i = 0
            while i < self.nbalt and self.alts_down[i] >= alt:
                i += 1
            if i == self.nbalt:  # to handle alt min
                i = self.nbalt - 1
            low_index = i  # grid low index
            high_index = i - 1  # grid high index
        return (high_index, low_index)

    def direct_loc_grid_h(self, row0, col0, steprow, stepcol, nbrow, nbcol, alt):
        """
        direct localization  grid at constant altitude
        TODO not tested.

        :param row0: grid origin (row)
        :type row0: int
        :param col0: grid origin (col)
        :type col0: int
        :param steprow: grid step (row)
        :type steprow: int
        :param stepcol: grid step (col)
        :type stepcol: int
        :param nbrow: grid nb row
        :type nbrow: int
        :param nbcol: grid nb col
        :type nbcol: int
        :param alt: altitude of the grid
        :type alt: float
        :return: direct localization grid
        :rtype: numpy.array
        """
        if isinstance(alt, (list, np.ndarray)):
            logging.warning("grid doesn't handle alt as array, first value is used")
            alt = alt[0]
        gldalt = np.zeros((3, nbrow, nbcol))
        (grid_index_up, grid_index_down) = self.return_grid_index(alt)
        alt_down = self.alts_down[grid_index_down]
        alt_up = self.alts_down[grid_index_up]
        alti_coef = (alt - alt_down) / (alt_up - alt_down)
        mats = [
            self.lon_data[grid_index_up : grid_index_down + 1, :, :],
            self.lat_data[grid_index_up : grid_index_down + 1, :, :],
        ]
        position = np.zeros(3)
        position[2] = alt
        for row_index in range(nbrow):
            row = row0 + steprow * row_index
            pos_row = (row - self.row0) / self.steprow
            for col_index in range(nbcol):
                col = col0 + stepcol * col_index
                pos_col = (col - self.col0) / self.stepcol
                # pylint disable for code clarity interpol_bilin_vectorized returns one list of 2 elements in this case
                # pylint: disable=unbalanced-tuple-unpacking
                [vlon, vlat] = interpol_bilin_grid(mats, self.nbrow, self.nbcol, pos_row, pos_col)
                position[0] = alti_coef * vlon[0] + (1 - alti_coef) * vlon[1]
                position[1] = alti_coef * vlat[0] + (1 - alti_coef) * vlat[1]
                gldalt[:, row_index, col_index] = position
        return gldalt

    # gitlab issue #58
    # pylint: disable=too-many-locals
    def estimate_inverse_loc_predictor(self, nbrow_pred=3, nbcol_pred=3):
        """
        initialize inverse localization polynomial predictor
        it composed of 4 polynoms estimated on a grid at hmin and hmax

        col_min = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        row_min = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        col_max = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        row_max = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        least squarred method is used to calculate coefficients, which are noramlized in [-1,1]

        :param nbrow_pred: predictor nb row (3 by default)
        :type nbrow_pred: int
        :param nbcol_pred: predictor nb col (3 by default)
        :type nbcol_pred: int
        """
        nb_alt = 2
        nb_coeff = 6
        nb_mes = nbcol_pred * nbrow_pred

        col_norm = np.linspace(-1.0, 1.0, nbcol_pred)
        row_norm = np.linspace(-1.0, 1.0, nbrow_pred)
        gcol_norm, grow_norm = np.meshgrid(col_norm, row_norm)
        glon, glat = self.interpolate_grid_in_altitude(nbrow_pred, nbcol_pred, nb_alt)

        # normalisation des variables
        (glon_min, glon_max) = (glon.min(), glon.max())
        (glat_min, glat_max) = (glat.min(), glat.max())
        lon_ofset = (glon_max + glon_min) / 2.0
        lon_scale = (glon_max - glon_min) / 2.0
        lat_ofset = (glat_max + glat_min) / 2.0
        lat_scale = (glat_max - glat_min) / 2.0

        glon_norm = (glon - lon_ofset) / lon_scale
        glat_norm = (glat - lat_ofset) / lat_scale

        col_ofset = (self.colmax + self.col0) / 2.0
        col_scale = (self.colmax - self.col0) / 2.0
        row_ofset = (self.rowmax + self.row0) / 2.0
        row_scale = (self.rowmax - self.row0) / 2.0

        glon2 = glon_norm * glon_norm
        glonlat = glon_norm * glat_norm
        glat2 = glat_norm * glat_norm

        a_min = np.zeros((nb_mes, nb_coeff))
        a_max = np.zeros((nb_mes, nb_coeff))
        b_col = np.zeros((nb_mes, 1))
        b_row = np.zeros((nb_mes, 1))

        # resolution des moindres carres
        imes = 0
        for irow in range(nbrow_pred):
            for icol in range(nbcol_pred):
                b_col[imes] = gcol_norm[irow, icol]
                b_row[imes] = grow_norm[irow, icol]
                a_min[imes, 0] = 1.0
                a_max[imes, 0] = 1.0
                a_min[imes, 1] = glon_norm[1, irow, icol]
                a_max[imes, 1] = glon_norm[0, irow, icol]
                a_min[imes, 2] = glat_norm[1, irow, icol]
                a_max[imes, 2] = glat_norm[0, irow, icol]
                a_min[imes, 3] = glon2[1, irow, icol]
                a_max[imes, 3] = glon2[0, irow, icol]
                a_min[imes, 4] = glat2[1, irow, icol]
                a_max[imes, 4] = glat2[0, irow, icol]
                a_min[imes, 5] = glonlat[1, irow, icol]
                a_max[imes, 5] = glonlat[0, irow, icol]
                imes += 1

        # Compute coefficients
        mat_a_min = np.array(a_min)
        mat_a_max = np.array(a_max)

        t_aa_min = mat_a_min.T @ mat_a_min
        t_aa_max = mat_a_max.T @ mat_a_max
        t_aa_min_inv = np.linalg.inv(t_aa_min)
        t_aa_max_inv = np.linalg.inv(t_aa_max)

        coef_col_min = t_aa_min_inv @ mat_a_min.T @ b_col
        coef_row_min = t_aa_min_inv @ mat_a_min.T @ b_row
        coef_col_max = t_aa_max_inv @ mat_a_max.T @ b_col
        coef_row_max = t_aa_max_inv @ mat_a_max.T @ b_row

        # seems not clear to understand with inverse_loc_predictor function ...
        self.pred_col_min = coef_col_min.flatten()
        self.pred_row_min = coef_row_min.flatten()
        self.pred_col_max = coef_col_max.flatten()
        self.pred_row_max = coef_row_max.flatten()
        self.pred_ofset_scale_lon = [lon_ofset, lon_scale]
        self.pred_ofset_scale_lat = [lat_ofset, lat_scale]
        self.pred_ofset_scale_row = [row_ofset, row_scale]
        self.pred_ofset_scale_col = [col_ofset, col_scale]

    def inverse_loc_predictor(self, lon, lat, alt=0.0):
        """
        evaluate inverse localization predictor at a given geographic position

        :param lon: longitude
        :type lon: float
        :param lat: latitude
        :type lat: float
        :param alt: altitude (0.0 by default)
        :type alt: float
        :return: sensor position and extrapolation state (row,col, is extrapolated)
        :rtype: tuple (float, float, boolean)
        """
        extrapolation_threshold = 20.0
        is_extrapolated = False
        altmin = self.alts_down[-1]
        altmax = self.alts_down[0]

        # normalization
        lon_n = (lon - self.pred_ofset_scale_lon[0]) / self.pred_ofset_scale_lon[1]
        lat_n = (lat - self.pred_ofset_scale_lat[0]) / self.pred_ofset_scale_lat[1]
        if abs(lon_n) > (1 + extrapolation_threshold / 100.0):
            logging.debug("Be careful: longitude extrapolation: %1.8f", lon_n)
            is_extrapolated = True
        if abs(lat_n) > (1 + extrapolation_threshold / 100.0):
            logging.debug("Be careful: latitude extrapolation: %1.8f", lat_n)
            is_extrapolated = True

        # polynome application
        vect_sol = np.array([1, lon_n, lat_n, lon_n**2, lat_n**2, lon_n * lat_n])
        col_min = ((self.pred_col_min * vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        row_min = ((self.pred_row_min * vect_sol).sum() * self.pred_ofset_scale_row[1]) + self.pred_ofset_scale_row[0]
        col_max = ((self.pred_col_max * vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        row_max = ((self.pred_row_max * vect_sol).sum() * self.pred_ofset_scale_row[1]) + self.pred_ofset_scale_row[0]

        h_x = (alt - altmin) / (altmax - altmin)
        col = (1 - h_x) * col_min + h_x * col_max
        row = (1 - h_x) * row_min + h_x * row_max

        return row, col, is_extrapolated

    # gitlab issue #58
    # pylint: disable=too-many-locals
    def inverse_partial_derivative(self, row, col, alt=0):
        """
        calculate partial derivative at a given geographic position
        it gives the matrix to apply to get sensor shifts from geographic ones
        it returns M matrix :
        [dcol,drow]T = M x [dlon,dlat]T
        dlon/dlat in microrad
        M is calculated on each node of grid
        M is necessary for direct localization inversion in iterative inverse loc

        :param lon: longitude
        :type lon: float
        :param lat: latitude
        :type lat: float
        :param alt: altitude (0.0 by default)
        :type alt: float
        :return: matrix
        :rtype: numpy.array
        """
        pos_row = (row - self.row0) / self.steprow
        pos_col = (col - self.col0) / self.stepcol
        index_row = np.floor(pos_row)
        index_col = np.floor(pos_col)

        index_row = min(index_row, self.nbrow - 2)
        index_row = max(index_row, 0)
        index_col = min(index_col, self.nbcol - 2)
        index_col = max(index_col, 0)
        index_row = int(index_row)
        index_col = int(index_col)

        (grid_index_up, grid_index_down) = self.return_grid_index(alt)

        lon_h00 = self.lon_data[grid_index_up, index_row, index_col]
        lon_h01 = self.lon_data[grid_index_up, index_row, index_col + 1]
        lon_h10 = self.lon_data[grid_index_up, index_row + 1, index_col]

        lon_b00 = self.lon_data[grid_index_down, index_row, index_col]
        lon_b01 = self.lon_data[grid_index_down, index_row, index_col + 1]
        lon_b10 = self.lon_data[grid_index_down, index_row + 1, index_col]

        lat_h00 = self.lat_data[grid_index_up, index_row, index_col]
        lat_h01 = self.lat_data[grid_index_up, index_row, index_col + 1]
        lat_h10 = self.lat_data[grid_index_up, index_row + 1, index_col]

        lat_b00 = self.lat_data[grid_index_down, index_row, index_col]
        lat_b01 = self.lat_data[grid_index_down, index_row, index_col + 1]
        lat_b10 = self.lat_data[grid_index_down, index_row + 1, index_col]

        dlon_ch = np.deg2rad(lon_h01 - lon_h00) / self.stepcol
        dlon_cb = np.deg2rad(lon_b01 - lon_b00) / self.stepcol
        dlon_lh = np.deg2rad(lon_h10 - lon_h00) / self.steprow
        dlon_lb = np.deg2rad(lon_b10 - lon_b00) / self.steprow

        dlat_ch = np.deg2rad(lat_h01 - lat_h00) / self.stepcol
        dlat_cb = np.deg2rad(lat_b01 - lat_b00) / self.stepcol
        dlat_lh = np.deg2rad(lat_h10 - lat_h00) / self.steprow
        dlat_lb = np.deg2rad(lat_b10 - lat_b00) / self.steprow

        h_x = (alt - self.alts_down[grid_index_down]) / (
            self.alts_down[grid_index_up] - self.alts_down[grid_index_down]
        )

        dlon_c = ((1 - h_x) * dlon_cb + (h_x) * dlon_ch) * 1e6
        dlat_c = ((1 - h_x) * dlat_cb + (h_x) * dlat_ch) * 1e6
        dlon_l = ((1 - h_x) * dlon_lb + (h_x) * dlon_lh) * 1e6
        dlat_l = ((1 - h_x) * dlat_lb + (h_x) * dlat_lh) * 1e6
        det = dlon_c * dlat_l - dlon_l * dlat_c
        partial_derivative_mat = np.array([[dlat_l, -dlon_l], [-dlat_c, dlon_c]]) / det
        if abs(det) <= 0.000000000001:
            logging.warning("inverse_loc() inverse_partial_derivative() determinant is null")
        return partial_derivative_mat

    def inverse_loc(self, lon, lat, alt=0.0, nb_iterations=15):
        """
        Inverse localization at a given geographic position
        First initialize position,
        * apply inverse predictor lon,lat,at ->  col_0,row_0
        * direct loc col_0,row_0 -> lon_0, lat_0
        Then iterative process:
        * calculate geographic error dlon,dlat
        * calculate sensor correction dlon,dlat -> dcol,drow
        * apply direct localization  -> lon_i,lat_i
        * compute localisation error
        * compute local inverse gradient
        * move along derivatives to compute row_i,col_i
        * loop until measurement  error is below threshold or max number of iterations


        TODO optimization (for loop,...)

        :param lon: longitude
        :type lon: float or 1D numpy.ndarray dtype=float64
        :param lat: latitude
        :type lat: float or 1D numpy.ndarray dtype=float64
        :param alt: altitude
        :type alt: float or 1D numpy.ndarray dtype=float64
        :param nb_iterations: max number of iterations (15 by default)
        :type nb_iterations: int
        :return: sensor position (row,col,alt)
        :rtype: tuple(1D np.array row position, 1D np.array col position, 1D np.array alt)
        """

        # Test added for rectification to work
        # TODO: refactoring with interfaces clean
        if not isinstance(lon, np.ndarray):
            lon = np.array([lon])
            lat = np.array([lat])
        if not isinstance(alt, np.ndarray):
            alt = np.array([alt])

        if alt.shape[0] != lon.shape[0]:
            alt = np.full(lon.shape[0], fill_value=alt[0])

        points_nb = len(lon)
        filter_nan = np.logical_not(np.logical_or(np.isnan(lon), np.isnan(lat)))
        points_nb_valid = np.sum(filter_nan)
        logging.debug("number of valid points %d,  number of points %d", points_nb_valid, points_nb)
        row = np.nan * np.zeros((points_nb))
        col = np.nan * np.zeros((points_nb))
        points_valid = np.arange(points_nb)[filter_nan]
        for point_index in points_valid:
            lon_i = lon[point_index]
            lat_i = lat[point_index]
            alt_i = alt[point_index]
            deg2mrad = np.deg2rad(1.0) * 1e6
            iteration = 0
            coslon = np.cos(np.deg2rad(lat_i))
            rtx = 1e-12 * 6378000**2
            row_i, col_i, _ = self.inverse_loc_predictor(lon_i, lat_i, alt_i)
            m2_error = 10.0
            valid_point = 0

            # Iterative process
            # while error in m2 > 1mm
            while (m2_error > 1e-6) and (iteration < nb_iterations):
                position = self.direct_loc_h(row_i, col_i, alt_i)
                dlon_microrad = (position[0][0] - lon_i) * deg2mrad
                dlat_microrad = (position[0][1] - lat_i) * deg2mrad
                m2_error = rtx * (dlat_microrad**2 + (dlon_microrad * coslon) ** 2)
                dsol = np.array([dlon_microrad, dlat_microrad])
                mat_dp = self.inverse_partial_derivative(row_i, col_i, alt_i)
                dimg = mat_dp @ dsol
                col_i += -dimg[0]
                row_i += -dimg[1]
                iteration += 1
                valid_point = 1
            if valid_point == 0:
                (row_i, col_i) = (None, None)
            row[point_index] = row_i
            col[point_index] = col_i
        return row, col, alt


def coloc(multi_h_grid_src, multi_h_grid_dst, dtm, origin, step, size):
    """
    colocalization grid on dtm
    localization on dtm from src grid, then inverse localization in right grid

    :param multi_h_grid_src: source grid
    :type multi_h_grid_src: shareloc.grid
    :param multi_h_grid_dst: destination grid
    :type multi_h_grid_dst: shareloc.grid
    :param origin: grid origin in src grid (row,col)
    :type origin: list(int)
    :param step: grid step (row,col)
    :type step: list(int)
    :param size: grid nb row and nb col
    :type size: list(int)
    :return: colocalization grid
    :rtype: numpy.array
    """
    [l0_src, c0_src] = origin
    [steprow_src, stepcol_src] = step
    [nbrow_src, nbcol_src] = size
    gricoloc = np.zeros((3, nbrow_src, nbcol_src))
    for index_row in range(nbrow_src):
        row = l0_src + steprow_src * index_row
        for index_col in range(nbcol_src):
            col = c0_src + stepcol_src * index_col
            # TODO: refacto interfaces to avoid np.squeeze if not necessary
            (lon, lat, alt) = np.squeeze(multi_h_grid_src.direct_loc_dtm(row, col, dtm))
            pos_dst = multi_h_grid_dst.inverse_loc(lon, lat, alt)
            pos_dst = np.array([*pos_dst])[:, 0]
            gricoloc[:, index_row, index_col] = pos_dst
    return gricoloc
