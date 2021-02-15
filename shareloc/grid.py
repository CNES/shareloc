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
localisations function from mutji h direct grids.
"""

import numpy as np
from shareloc.readwrite import read_bsq_hd
from shareloc.math_utils import interpol_bilin, interpol_bilin_vectorized


class Grid:
    """multi H direct localization grid handling class"""

    # gitlab issue #58
    # pylint: disable=too-many-instance-attributes
    def __init__(self, grid_filename, grid_format="bsq"):
        """
        Constructor
        :param grid_filename: grid filename
        :type grid_filename: string
        :param grid_format: grid format (by default bsq)
        :type grid_format: string
        """
        self.filename = grid_filename
        self.format = grid_format
        self.row0 = None
        self.col0 = None
        self.nbrow = None
        self.nbcol = None
        self.steprow = None
        self.stepcol = None
        self.repter = None
        self.nbalt = None
        self.index_alt = {}
        self.gld_lon = None
        self.gld_lat = None
        self.alts_down = []
        self.rowmax = None
        self.colmax = None
        self.load()
        self.type = "multi H grid"

    def load(self):
        """
        header and grid loading function
        2 data cubes are defined :
        - gld_lon : [alt,row,col]
        - gld_lat : [alt,row,col]
        bsq grids are stored by increasing altitude H0 ... Hx
        internal structure is decreasing one
        """
        if self.format == "bsq":
            dico_a_lire = {
                "nbrow": ("LINES", int),
                "nbcol": ("COLUMNS", int),
                "bpp": ("BITS PER PIXEL", int),
                "nbalt": ("NB ALT", int),
                "stepcol": ("PAS COL", float),
                "steprow": ("PAS LIG", float),
                "col0": ("COL0", float),
                "row0": ("LIG0", float),
                "repter": ("REFERENTIEL TERRESTRE", str),
            }

            nom_hd = self.filename[:-4] + "1.hd"
            dico_hd = read_bsq_hd(nom_hd, dico_a_lire)

            for var in dico_hd:
                setattr(self, var, dico_hd[var])

            # """renvoie une structure 3D [i_alt][l,c]"""
            gld_lon = np.zeros((self.nbalt, self.nbrow, self.nbcol))
            gld_lat = np.zeros((self.nbalt, self.nbrow, self.nbcol))

            codage = float

            for alt_layer in range(self.nbalt):
                inverse_index = self.nbalt - alt_layer
                nom_gri_lon = self.filename[:-4] + str(inverse_index) + ".c1"
                nom_gri_lat = self.filename[:-4] + str(inverse_index) + ".c2"
                nom_hd = self.filename[:-4] + str(inverse_index) + ".hd"

                gld_lon[alt_layer, :, :] = np.fromfile(nom_gri_lon, dtype=codage).reshape((self.nbrow, self.nbcol))
                gld_lat[alt_layer, :, :] = np.fromfile(nom_gri_lat, dtype=codage).reshape((self.nbrow, self.nbcol))

                dico_hd = read_bsq_hd(nom_hd, {"index": ("ALT INDEX", int), "alt": ("ALTITUDE", float)})
                self.index_alt[dico_hd["index"]] = dico_hd["alt"]

            self.gld_lon = gld_lon
            self.gld_lat = gld_lat
            self.alts_down = [self.index_alt[_] for _ in range(int(self.nbalt - 1), -1, -1)]
            self.rowmax = self.row0 + self.steprow * (self.nbrow - 1)
            self.colmax = self.col0 + self.stepcol * (self.nbcol - 1)
        else:
            print("dtm format is not handled")

    def get_alt_min_max(self):
        """
        returns altitudes min and max layers
        :return alt_min,lat_max
        :rtype list
        """
        return [self.alts_down[-1], self.alts_down[0]]

    def direct_loc_h(self, row, col, alt, fill_nan=False):
        """
        direct localization at constant altitude
        :param row :  line sensor position
        :type row : float or 1D numpy.ndarray dtype=float64
        :param col :  column sensor position
        :type col : float or 1D numpy.ndarray dtype=float64
        :param alt :  altitude
        :type alt : float
        :param fill_nan : fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan : boolean
        :return ground position (lon,lat,h)
        :rtype numpy.ndarray
        """
        if fill_nan:
            print("fill nan {}".format(fill_nan))
        # faire une controle sur row / col !!!!
        # 0.5 < row < rowmax
        (grid_index_up, grid_index_down) = self.return_grid_index(alt)
        alt_down = self.alts_down[grid_index_down]
        alt_up = self.alts_down[grid_index_up]
        alti_coef = (alt - alt_down) / (alt_up - alt_down)
        mats = [
            self.gld_lon[grid_index_up : grid_index_down + 1, :, :],
            self.gld_lat[grid_index_up : grid_index_down + 1, :, :],
        ]

        if not isinstance(col, (list, np.ndarray)):
            col = np.array([col])
            row = np.array([row])

        position = np.zeros((col.size, 3))
        position[:, 2] = alt
        pos_row = (row - self.row0) / self.steprow
        pos_col = (col - self.col0) / self.stepcol
        # pylint disable for code clarity interpol_bilin_vectorized returns one list of 2 elements in this case
        # pylint: disable=unbalanced-tuple-unpacking
        [vlon, vlat] = interpol_bilin_vectorized(mats, self.nbrow, self.nbcol, pos_row, pos_col)
        position[:, 0] = alti_coef * vlon[0, :] + (1 - alti_coef) * vlon[1, :]
        position[:, 1] = alti_coef * vlat[0, :] + (1 - alti_coef) * vlat[1, :]

        return np.squeeze(position)

    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :param dtm : dtm model
        :type dtm  : shareloc.dtm
        :return ground position (lon,lat,h)
        :rtype numpy.array
        """
        los = np.zeros((3, self.nbalt))
        loslonlat = self.interpolate_grid_in_plani(row, col)
        los[0, :] = loslonlat[0]
        los[1, :] = loslonlat[1]
        los[2, :] = self.alts_down
        los = los.T
        (__, __, point_b, alti) = dtm.intersect_dtm_cube(los)
        (__, __, point_dtm) = dtm.intersection(los, point_b, alti)
        return point_dtm

    def los_extrema(self, row, col, alt_min, alt_max):
        """
        compute los extrema
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :param alt_min : los alt min
        :type alt_min  : float
        :param alt_max : los alt max
        :type alt_max : float
        :return los extrema
        :rtype numpy.array (2x3)
        """
        los_edges = np.zeros([2, 3])
        los_edges[0, :] = self.direct_loc_h(row, col, alt_max)
        los_edges[1, :] = self.direct_loc_h(row, col, alt_min)
        return los_edges

    def interpolate_grid_in_plani(self, row, col):
        """
        interpolate positions on multi h grid
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :return interpolated positions
        :rtype list
        """
        pos_row = (row - self.row0) / self.steprow
        pos_col = (col - self.col0) / self.stepcol
        mats = [self.gld_lon, self.gld_lat]
        res = interpol_bilin(mats, self.nbrow, self.nbcol, pos_row, pos_col)
        return res

    def interpolate_grid_in_altitude(self, nbrow, nbcol, nbalt=None):
        """
        interpolate equally spaced grid (in altitude)
        :param nbrow :  grid nb row
        :type nbrow : int
        :param nbcol :  grid nb col
        :type nbcol : int
        :param nbalt :  grid nb alt, of None self.nbalt is used instead
        :type nbalt : int
        :return equally spaced grid
        :rtype numpy.array
        """
        if not nbalt:
            nbalt = self.nbalt
            list_alts = self.alts_down
        else:
            list_alts = np.linspace(self.alts_down[0], self.alts_down[-1], nbalt)

        gld_lon = np.zeros((nbalt, nbrow, nbcol))
        gld_lat = np.zeros((nbalt, nbrow, nbcol))
        # """genere un cube de visee interpole de nrow/ncol visee"""
        # row_max = self.row0 + self.steprow * (self.nbrow-1)
        # col_max = self.col0 + self.stepcol * (self.nbcol-1)

        steprow = (self.rowmax - self.row0) / (nbrow - 1)
        stepcol = (self.colmax - self.col0) / (nbcol - 1)

        for index, alt in enumerate(list_alts):

            res = self.direct_loc_grid_h(self.row0, self.col0, steprow, stepcol, nbrow, nbcol, alt)
            gld_lon[index] = res[0]
            gld_lat[index] = res[1]
        return gld_lon, gld_lat

    def direct_loc_grid_dtm(self, row0, col0, steprow, stepcol, nbrow, nbcol, dtm):
        """
        direct localization  grid on dtm
        :param row0 :  grid origin (row)
        :type row0 : int
        :param col0 :  grid origin (col)
        :type col0 : int
        :param steprow :  grid step (row)
        :type steprow : int
        :param stepcol :  grid step (col)
        :type stepcol : int
        :param nbrow :  grid nb row
        :type nbrow : int
        :param nbcol :  grid nb col
        :type nbcol : int
        :param dtm : dtm model
        :type dtm  : shareloc.dtm
        :return direct localization grid
        :rtype numpy.array
        """
        glddtm = np.zeros((3, nbrow, nbcol))
        los = np.zeros((3, self.nbalt))
        for i in range(nbrow):
            for j in range(nbcol):
                col = col0 + stepcol * j
                row = row0 + steprow * i
                loslonlat = self.interpolate_grid_in_plani(row, col)
                los[0, :] = loslonlat[0]
                los[1, :] = loslonlat[1]
                los[2, :] = self.alts_down
                los_t = los.T
                (__, __, point_b, alti) = dtm.intersect_dtm_cube(los_t)
                (__, __, point_r) = dtm.intersection(los_t, point_b, alti)
                glddtm[:, i, j] = point_r
        return glddtm

    def return_grid_index(self, alt):
        """
        return layer index enclosing a given altitude
        :param alt :  altitude
        :type alt : float
        :return grid index (up,down)
        :rtype tuple
        """
        if alt > self.alts_down[0]:
            (indicehaut, indicebas) = (0, 0)
        elif alt < self.alts_down[-1]:
            (indicehaut, indicebas) = (self.nbalt - 1, self.nbalt - 1)
        else:
            i = 0
            while i < self.nbalt and self.alts_down[i] >= alt:
                i += 1
            if i == self.nbalt:  # pour gerer alt min
                i = self.nbalt - 1
            indicebas = i  # indice grille bas
            indicehaut = i - 1  # indice grille haut
        return (indicehaut, indicebas)

    def direct_loc_grid_h(self, row0, col0, steprow, stepcol, nbrow, nbcol, alt):
        """
        direct localization  grid at constant altitude
        :param row0 :  grid origin (row)
        :type row0 : int
        :param col0 :  grid origin (col)
        :type col0 : int
        :param steprow :  grid step (row)
        :type steprow : int
        :param stepcol :  grid step (col)
        :type stepcol : int
        :param nbrow :  grid nb row
        :type nbrow : int
        :param nbcol :  grid nb col
        :type nbcol : int
        :param alt : altitude of the grid
        :type alt  : float
        :return direct localization grid
        :rtype numpy.array
        """
        gldalt = np.zeros((3, nbrow, nbcol))
        (grid_index_up, grid_index_down) = self.return_grid_index(alt)
        alt_down = self.alts_down[grid_index_down]
        alt_up = self.alts_down[grid_index_up]
        alti_coef = (alt - alt_down) / (alt_up - alt_down)
        mats = [
            self.gld_lon[grid_index_up : grid_index_down + 1, :, :],
            self.gld_lat[grid_index_up : grid_index_down + 1, :, :],
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
                [vlon, vlat] = interpol_bilin(mats, self.nbrow, self.nbcol, pos_row, pos_col)
                position[0] = alti_coef * vlon[0] + (1 - alti_coef) * vlon[1]
                position[1] = alti_coef * vlat[0] + (1 - alti_coef) * vlat[1]
                gldalt[:, row_index, col_index] = position
        return gldalt

    # gitlab issue #58
    # pylint: disable=too-many-locals
    def estimate_inverse_loc_predictor(self, nbrow_pred=3, nbcol_pred=3):
        """
        initialize inverse localization polynomial predictor
        it composed of 4 polynoms estimated on a grid at hmin and hmax :

        col_min = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        row_min = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        col_max = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        row_max = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        least squarred method is used to calculate coefficients, which are noramlized in [-1,1]
        :param nbrow_pred :  predictor nb row (3 by default)
        :type nbrow_pred : int
        :param nbcol_pred :  predictor nb col (3 by default)
        :type nbcol_pred : int
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

        # Calcul des coeffcients
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

        setattr(self, "pred_col_min", coef_col_min.flatten())
        setattr(self, "pred_row_min", coef_row_min.flatten())
        setattr(self, "pred_col_max", coef_col_max.flatten())
        setattr(self, "pred_row_max", coef_row_max.flatten())
        setattr(self, "pred_ofset_scale_lon", [lon_ofset, lon_scale])
        setattr(self, "pred_ofset_scale_lat", [lat_ofset, lat_scale])
        setattr(self, "pred_ofset_scale_row", [row_ofset, row_scale])
        setattr(self, "pred_ofset_scale_col", [col_ofset, col_scale])

    def inverse_loc_predictor(self, lon, lat, alt=0.0):
        """
        evaluate inverse localization predictor at a given geographic position
        :param lon : longitude
        :type lon : float
        :param lat : latitude
        :type lat : float
        :param alt : altitude (0.0 by default)
        :type alt : float
        :return sensor position (row,col, is extrapolated)
        :rtype tuple (float,float,boolean)
        """
        seuil_extrapol = 20.0
        extrapol = False
        altmin = self.alts_down[-1]
        altmax = self.alts_down[0]

        # normalisation
        lon_n = (lon - self.pred_ofset_scale_lon[0]) / self.pred_ofset_scale_lon[1]
        lat_n = (lat - self.pred_ofset_scale_lat[0]) / self.pred_ofset_scale_lat[1]
        if abs(lon_n) > (1 + seuil_extrapol / 100.0):
            # print "Attention, en extrapolation de 20% en longitude:",lon_n
            extrapol = True
        if abs(lat_n) > (1 + seuil_extrapol / 100.0):
            # print "Attention, en extrapolation de 20% en latitude:",lat_n
            extrapol = True

        # application polynome
        vect_sol = np.array([1, lon_n, lat_n, lon_n ** 2, lat_n ** 2, lon_n * lat_n])
        col_min = ((self.pred_col_min * vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        row_min = ((self.pred_row_min * vect_sol).sum() * self.pred_ofset_scale_row[1]) + self.pred_ofset_scale_row[0]
        col_max = ((self.pred_col_max * vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        row_max = ((self.pred_row_max * vect_sol).sum() * self.pred_ofset_scale_row[1]) + self.pred_ofset_scale_row[0]

        h_x = (alt - altmin) / (altmax - altmin)
        col = (1 - h_x) * col_min + h_x * col_max
        row = (1 - h_x) * row_min + h_x * row_max

        if row > self.rowmax:
            row = self.rowmax
        if row < self.row0:
            row = self.row0
        if col > self.colmax:
            col = self.colmax
        if col < self.col0:
            col = self.col0
        return row, col, extrapol

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
        :param lon : longitude
        :type lon : float
        :param lat : latitude
        :type lat : float
        :param alt : altitude (0.0 by default)
        :type alt : float
        :return matrix
        :rtype numpy.array
        """
        pos_row = (row - self.row0) / self.steprow
        pos_col = (col - self.col0) / self.stepcol
        index_row = int(np.floor(pos_row))
        index_col = int(np.floor(pos_col))
        (grid_index_up, grid_index_down) = self.return_grid_index(alt)

        lon_h00 = self.gld_lon[grid_index_up, index_row, index_col]
        lon_h01 = self.gld_lon[grid_index_up, index_row, index_col + 1]
        lon_h10 = self.gld_lon[grid_index_up, index_row + 1, index_col]

        lon_b00 = self.gld_lon[grid_index_down, index_row, index_col]
        lon_b01 = self.gld_lon[grid_index_down, index_row, index_col + 1]
        lon_b10 = self.gld_lon[grid_index_down, index_row + 1, index_col]

        lat_h00 = self.gld_lat[grid_index_up, index_row, index_col]
        lat_h01 = self.gld_lat[grid_index_up, index_row, index_col + 1]
        lat_h10 = self.gld_lat[grid_index_up, index_row + 1, index_col]

        lat_b00 = self.gld_lat[grid_index_down, index_row, index_col]
        lat_b01 = self.gld_lat[grid_index_down, index_row, index_col + 1]
        lat_b10 = self.gld_lat[grid_index_down, index_row + 1, index_col]

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
        if abs(det) > 0.000000000001:
            partial_derivative_mat = np.array([[dlat_l, -dlon_l], [-dlat_c, dlon_c]]) / det
        else:
            print("nul determinant")
        return partial_derivative_mat

    def inverse_loc(self, lon, lat, alt=0.0, nb_iterations=15):
        """
        inverse localization at a given geographic position
        First initialize position,
        * apply inverse predictor lon,lat,at ->  col_0,row_0
        * direct loc col_0,row_0 -> lon_0, lat_0
        Then iterative process:
        * calculate geographic error dlon,dlat
        * calculate senor correction dlon,dlat -> dcol,drow
        * apply direct localization  -> lon_i,lat_i
        :param lon : longitude
        :type lon: float
        :param lat : latitude
        :type lat: float
        :param alt : altitude
        :type alt: float
        :param nb_iterations : max number of iterations (15 by default)
        :type nb_iterations : int
        :return sensor position (row,col,alt, is_valid)
        :rtype tuple (float,float,float,boolean)
        """

        deg2mrad = np.deg2rad(1.0) * 1e6
        iteration = 0
        coslon = np.cos(np.deg2rad(lat))
        rtx = 1e-12 * 6378000 ** 2
        (row_i, col_i, extrapol) = self.inverse_loc_predictor(lon, lat, alt)
        erreur_m2 = 10.0
        point_valide = 0
        if not extrapol:
            # Processus iteratif
            # while erreur > seuil:1mm
            while (erreur_m2 > 1e-6) and (iteration < nb_iterations):
                # print k,row_i,col_i
                position = self.direct_loc_h(row_i, col_i, alt)
                dlon_microrad = (position[0] - lon) * deg2mrad
                dlat_microrad = (position[1] - lat) * deg2mrad
                erreur_m2 = rtx * (dlat_microrad ** 2 + (dlon_microrad * coslon) ** 2)
                dsol = np.array([dlon_microrad, dlat_microrad])
                mat_dp = self.inverse_partial_derivative(row_i, col_i, alt)
                dimg = mat_dp @ dsol
                col_i += -dimg[0]
                row_i += -dimg[1]
                iteration += 1
                point_valide = 1
        return (row_i, col_i, alt, point_valide)


def coloc(multi_h_grid_src, multi_h_grid_dst, dtm, origin, step, size):
    """
    colocalization grid on dtm
    localization on dtm from src grid, then inverse localization in right grid
    :param multi_h_grid_src : source grid
    :type multi_h_grid_src : shareloc.grid
    :param multi_h_grid_dst : destination grid
    :type multi_h_grid_dst : shareloc.grid
     :param origin :  grid origin in src grid (row,col)
     :type origin : list(int)
     :param step :  grid step (row,col)
     :type step : list(int)
     :param size :  grid nb row and nb col
     :type size : list(int)
     :return colocalization grid
     :rtype numpy.array
    """
    [l0_src, c0_src] = origin
    [steprow_src, stepcol_src] = step
    [nbrow_src, nbcol_src] = size
    gricoloc = np.zeros((4, nbrow_src, nbcol_src))
    for index_row in range(nbrow_src):
        row = l0_src + steprow_src * index_row
        for index_col in range(nbcol_src):
            col = c0_src + stepcol_src * index_col
            (lon, lat, alt) = multi_h_grid_src.direct_loc_dtm(row, col, dtm)
            pos_dst = multi_h_grid_dst.inverse_loc(lon, lat, alt)
            gricoloc[:, index_row, index_col] = pos_dst
    return gricoloc
