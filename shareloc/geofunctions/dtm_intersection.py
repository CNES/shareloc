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
This module contains the DTMIntersection class to handle DTM intersection.
"""

# Standard imports
import logging

# Third party imports
import numpy as np
from scipy import interpolate

# Shareloc imports
from shareloc.dtm_image import DTMImage
from shareloc.image import Image
from shareloc.math_utils import interpol_bilin
from shareloc.proj_utils import coordinates_conversion


def interpolate_geoid_height(geoid_filename, positions, interpolation_method="linear"):
    """
    terrain to index conversion
    retrieve geoid height above ellispoid

    :param geoid_filename: geoid_filename
    :type geoid_filename: str
    :param positions: geodetic coordinates
    :type positions: 2D numpy array: (number of points, [long coord, lat coord])
    :param interpolation_method: default is 'linear' (interpn interpolation method)
    :type interpolation_method: str
    :return: geoid height
    :rtype: 1 numpy array (number of points)
    """

    geoid_image = Image(geoid_filename, read_data=True)

    # Prepare grid for interpolation
    row_indexes = np.arange(0, geoid_image.nb_rows, 1)
    col_indexes = np.arange(0, geoid_image.nb_columns, 1)
    points = (row_indexes, col_indexes)

    # add modulo lon/lat
    min_lon = geoid_image.origin_col
    max_lon = min_lon + geoid_image.nb_columns * geoid_image.pixel_size_col
    positions[:, 0] += (positions[:, 0] + min_lon < 0) * 360.0
    positions[:, 0] -= (positions[:, 0] - max_lon > 0) * 360.0
    if np.any(np.abs(positions[:, 1]) > 90.0):
        raise RuntimeError("Geoid cannot handle latitudes greater than 90 deg.")
    indexes_geoid = geoid_image.transform_physical_point_to_index(positions[:, 1], positions[:, 0])
    return interpolate.interpn(points, geoid_image.data[:, :], indexes_geoid, method=interpolation_method)


class DTMIntersection:
    """
    DTMIntersection class dedicated to earth elevation handling with a DTM

    we work in cell convention [0,0] is the first cell center (not [0.5,0.5])
    """

    # gitlab issue #56
    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        dtm_filename,
        geoid_filename=None,
        roi=None,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    ):
        """
        Constructor
        :param dtm_filename: dtm filename
        :type dtm_filename: string
        :param geoid_filename: geoid filename, if None datum is ellispoid
        :type geoid_filename: string
        :param roi: region of interest [row_min,col_min,row_max,col_max] or [xmin,y_min,x_max,y_max] if
             roi_is_in_physical_space activated
        :type roi: list
        :param roi_is_in_physical_space: roi value in physical space
        :type roi_is_in_physical_space: bool
        :param fill_nodata:  fill_nodata strategy in None/'constant'/'min'/'median'/'max'/'mean'/'rio_fillnodata'/
        :type fill_nodata: str
        :param fill_value:  fill value for constant strategy. fill value is used for 'roi_fillnodata' residuals nodata,
            if None 'min' is used
        :type fill_value: float
        """
        self.dtm_file = dtm_filename
        self.alt_data = None
        self.alt_min = None
        self.alt_max = None
        self.origin_x = None
        self.origin_y = None
        self.pixel_size_x = None
        self.pixel_size_y = None
        self.plane_coef_a = None
        self.plane_coef_b = None
        self.plane_coef_c = None
        self.plane_coef_d = None
        self.alt_min_cell = None
        self.alt_max_cell = None
        self.tol_z = 0.0001

        # DTM reading
        datum = "ellipsoid"
        if geoid_filename is not None:
            datum = "geoid"

        self.dtm_image = DTMImage(
            self.dtm_file,
            read_data=True,
            roi=roi,
            roi_is_in_physical_space=roi_is_in_physical_space,
            datum=datum,
            fill_nodata=fill_nodata,
            fill_value=fill_value,
        )
        self.epsg = self.dtm_image.epsg
        self.alt_data = self.dtm_image.data[:, :].astype("float64")
        if self.dtm_image.datum == "geoid":
            logging.debug("remove geoid height")
            if geoid_filename is not None:
                self.grid_row, self.grid_col = np.mgrid[
                    0 : self.dtm_image.nb_rows : 1, 0 : self.dtm_image.nb_columns : 1
                ]
                lat, lon = self.dtm_image.transform_index_to_physical_point(self.grid_row, self.grid_col)
                positions = np.vstack([lon.flatten(), lat.flatten()]).transpose()
                if self.epsg != 4326:
                    positions = coordinates_conversion(positions, self.epsg, 4326)
                geoid_height = interpolate_geoid_height(geoid_filename, positions)
                self.alt_data += geoid_height.reshape(lon.shape)
            else:
                logging.warning("DTM datum is geoid but no geoid file is given")
        else:
            logging.debug("no geoid file is given dtm is assumed to be w.r.t ellipsoid")

        self.init_min_max()
        self.alt_max = self.alt_data.max()
        self.alt_min = self.alt_data.min()

        self.plane_coef_a = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.plane_coef_b = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        self.plane_coef_c = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        self.plane_coef_d = np.array(
            [0.0, self.dtm_image.nb_rows - 1.0, 0.0, self.dtm_image.nb_columns - 1.0, self.alt_min, self.alt_max]
        )

        self.plans = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, self.dtm_image.nb_rows - 1.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, self.dtm_image.nb_columns - 1.0],
                [0.0, 0.0, 1.0, self.alt_min],
                [0.0, 0.0, 1.0, self.alt_max],
            ]
        )

    def eq_plan(self, i, position):
        """
        return evaluation of equation on a plane on DTM cube

        :param i: face index
        :type i: int
        :param position: position
        :type position: numpy.array (1x3)
        :return: evaluation on the plan
        :rtype: float
        """
        return (
            self.plane_coef_a[i] * position[0]
            + self.plane_coef_b[i] * position[1]
            + self.plane_coef_c[i] * position[2]
            - self.plane_coef_d[i]
        )

    def ter_to_index(self, vect_ter):
        """
        terrain to index conversion

        :param vect_ter: terrain coordinate (lon,lat)
        :type vect_ter: array (1x2 or 1x3) if dimension is 3 , last coordinate is unchanged (alt)
        :return: index coordinates (col,row)
        :rtype: array (1x2 or 1x3)
        """
        vect_dtm = vect_ter.copy()
        (vect_dtm[0], vect_dtm[1]) = self.dtm_image.transform_physical_point_to_index(vect_ter[1], vect_ter[0])
        return vect_dtm

    def ters_to_indexs(self, vect_ters):
        """
        terrain to index conversion

        :param vect_ters: terrain coordinates
        :type vect_ters: array (nx2 or nx3) if dimension is 3 , last coordinates is unchanged (alt)
        :return: index coordinates
        :rtype: array (nx2 or nx3)
        """
        vect_dtms = vect_ters.copy()
        for i, vect_ter in enumerate(vect_ters):
            vect_dtms[i, :] = self.ter_to_index(vect_ter)
        return vect_dtms

    def index_to_ter(self, vect_dtm):
        """
        index to terrain conversion

        :param vect_dtm: index coordinate (col,row)
        :type vect_dtm: array (1x2 or 1x3) if dimension is 3 , last coordinate is unchanged (alt)
        :return: terrain coordinates (lon,lat)
        :rtype: array (1x2 or 1x3)
        """
        vect_ter = vect_dtm.copy()
        (vect_ter[1], vect_ter[0]) = self.dtm_image.transform_index_to_physical_point(vect_dtm[0], vect_dtm[1])
        return vect_ter

    def get_alt_offset(self, epsg):
        """
        returns min/amx altitude offset between dtm coordinates system and another one

        :param epsg: epsg code to compare with
        :type epsg: int
        :return: min/max altimetric difference between epsg in parameter minus dtm alti expressed in dtm epsg
        :rtype: list of float (1x2)
        """
        if epsg != self.epsg:
            alti_moy = (self.alt_min + self.alt_max) / 2.0
            corners = np.zeros([4, 2])
            corners[:, 0] = [0.0, 0.0, self.dtm_image.nb_rows, self.dtm_image.nb_rows]
            corners[:, 1] = [0.0, self.dtm_image.nb_columns, self.dtm_image.nb_columns, 0.0]
            corners -= 0.5  # index to corner
            ground_corners = np.zeros([3, 4])
            ground_corners[2, :] = alti_moy
            ground_corners[1::-1, :] = self.dtm_image.transform_index_to_physical_point(corners[:, 0], corners[:, 1])
            converted_corners = coordinates_conversion(ground_corners.transpose(), self.epsg, epsg)
            return [np.min(converted_corners[:, 2]) - alti_moy, np.max(converted_corners[:, 2]) - alti_moy]
        return [0.0, 0.0]

    def interpolate(self, pos_row, pos_col):
        """
        interpolate altitude
        if interpolation is done outside DTM the penultimate index is used (or first if d is negative).

        :param pos_row: cell position row
        :type pos_row: float
        :param pos_col: cell position col
        :type pos_col: float
        :return: interpolated altitude
        :rtype: float
        """
        alt = interpol_bilin(
            [self.alt_data[np.newaxis, :, :]], self.dtm_image.nb_rows, self.dtm_image.nb_columns, pos_row, pos_col
        )[0][0]
        return alt

    def init_min_max(self):
        """
        initialize min/max at each dtm cell
        """
        column_nb_minus = self.dtm_image.nb_columns - 1
        row_nb_minus = self.dtm_image.nb_rows - 1

        self.alt_max_cell = np.zeros((row_nb_minus, column_nb_minus))
        self.alt_min_cell = np.zeros((row_nb_minus, column_nb_minus))

        # Computation of the min and max altitudes of the DTM
        n_min = 32000
        n_max = -32000

        for i in range(row_nb_minus):
            k = 0
            for j in range(column_nb_minus):
                k += 1
                d_alt_min = n_min
                d_alt_max = n_max

                d_z1 = self.alt_data[i, j]
                if d_z1 < d_alt_min:
                    d_alt_min = d_z1
                if d_z1 > d_alt_max:
                    d_alt_max = d_z1

                d_z2 = self.alt_data[i, j + 1]
                if d_z2 < d_alt_min:
                    d_alt_min = d_z2
                if d_z2 > d_alt_max:
                    d_alt_max = d_z2

                d_z3 = self.alt_data[i + 1, j]
                if d_z3 < d_alt_min:
                    d_alt_min = d_z3
                if d_z3 > d_alt_max:
                    d_alt_max = d_z3

                d_z4 = self.alt_data[i + 1, j + 1]
                if d_z4 < d_alt_min:
                    d_alt_min = d_z4
                if d_z4 > d_alt_max:
                    d_alt_max = d_z4

                # GDN Correction BUG Intersector
                # It is important not to take the mathematical rounding for several reasons:
                # 1. the subsequent algorithm does not always resist well when the mesh is flat,
                #    it is therefore not recommended to output i_altmin = i_altmax
                #    unless this is really the case in real values
                # 2. the consecutive meshes must not be initialized in the same way by rounding
                #   because if the min altitude of one corresponds to
                # the max altitude of the other they must be distinguished
                #     by a ceil and a floor so that the cubes overlap slightly in altitude
                #    and not strictly contiguous
                i_altmin = np.floor(d_alt_min)
                i_altmax = np.ceil(d_alt_max)
                self.alt_min_cell[i, j] = i_altmin
                self.alt_max_cell[i, j] = i_altmax

    # gitlab issue #56
    # pylint: disable=too-many-branches
    def intersect_dtm_cube(self, los):  # noqa: C901
        """
        DTM cube intersection

        :param los:  line of sight
        :type los: numpy.array
        :return: intersection information
            (True,an intersection has been found ?, (lon,lat) of dtm position, altitude)
        :rtype: tuple (bool, bool, numpy.array, float)
        """
        # los: (n,3):
        point_b = None
        h_intersect = None
        (coord_col_i, coord_row_i, coord_alt_i, alti_layer_i) = ([], [], [], [])
        nbalt = los.shape[0]
        los_index = self.ters_to_indexs(los)
        # -----------------------------------------------------------------------
        # Number of valid intersections found
        nbi = 0
        # -----------------------------------------------------------------------
        # We loop on the 6 planes of the DTM cube
        for plane_index in range(6):
            # -----------------------------------------------------------------------
            # Init the vertex of the geometric line of sight
            los_hat = los_index[0, :]
            # -----------------------------------------------------------------------
            # Init position parallel to the plane
            # print self.plans[plane_index,:-1]
            # los_hat_onplane = (self.plans[plane_index,:-1]*los_hat).sum() - self.d[plane_index]
            # print los_hat_onplane
            los_hat_onplane = self.eq_plan(plane_index, los_hat)
            # -----------------------------------------------------------------------
            # Loop on line of sight segments
            # and control if we cross or not the current face plane_index of DTM cube
            for alti_layer in range(nbalt):
                # -----------------------------------------------------------------------
                # Transfer point B into point A
                los_a = los_hat_onplane
                s_a = los_hat.copy()
                # -----------------------------------------------------------------------
                # Reinit point B
                los_hat = los_index[alti_layer, :]  # we iterate on different los points
                # print los_hat
                # -----------------------------------------------------------------------
                # Init position parallel to the plane
                los_hat_onplane = self.eq_plan(plane_index, los_hat)
                # print 'posAposB',los_a,los_hat_onplane
                # -----------------------------------------------------------------------
                # Intersection test: los_a and los_hat_onplane with opposite sign
                if los_a * los_hat_onplane <= 0:
                    if not los_a and not los_hat_onplane:
                        # Too many solutions !! (A and B on the same plane) # How to deal with it ?
                        logging.warning("too many solutions in DTM intersection")
                    # -----------------------------------------------------------------------
                    if not los_a:
                        # A is solution (it is on the same plane)
                        coord_col_i.append(s_a[0])
                        coord_row_i.append(s_a[1])
                        coord_alt_i.append(s_a[2])
                        alti_layer_i.append(alti_layer - 1)
                    # -----------------------------------------------------------------------
                    elif not los_hat_onplane:
                        # B is solution (it is on the same plane)
                        coord_col_i.append(los_hat[0])
                        coord_row_i.append(los_hat[1])
                        coord_alt_i.append(los_hat[2])
                        alti_layer_i.append(alti_layer)
                    # -----------------------------------------------------------------------
                    else:
                        # -----------------------------------------------------------------------
                        # A and B are on either side of the plane
                        # Intersection interpolation coefficients between A and B
                        interp_coef_a = los_hat_onplane / (los_hat_onplane - los_a)
                        interp_coef_b = -los_a / (los_hat_onplane - los_a)
                        # Assignment or interpolation
                        # NB : to avoid test problems
                        #      <BeOnCube> (see further)
                        # . coordinate <u> (line)
                        # -----------------------------------------------------------------------
                        if plane_index < 2:
                            coord_col_i.append(self.plane_coef_d[plane_index])
                        # -----------------------------------------------------------------------
                        else:
                            coord_col_i.append(interp_coef_a * s_a[0] + interp_coef_b * los_hat[0])
                        # -----------------------------------------------------------------------
                        # . coordinate <v> (column)
                        if 1 < plane_index < 4:
                            coord_row_i.append(self.plane_coef_d[plane_index])
                        else:
                            coord_row_i.append(interp_coef_a * s_a[1] + interp_coef_b * los_hat[1])
                        # -----------------------------------------------------------------------
                        # . coordinate <z> (altitude)
                        if plane_index > 3:
                            coord_alt_i.append(self.plane_coef_d[plane_index])
                        # -----------------------------------------------------------------------
                        else:
                            coord_alt_i.append(interp_coef_a * s_a[2] + interp_coef_b * los_hat[2])
                        # . coordinate <h> (line of sight  x-axis)
                        alti_layer_i.append(alti_layer - interp_coef_a)  # index non entier de l'intersection
                    # -----------------------------------------------------------------------
                    # Incrementing the number of intersections found
                    nbi += 1
                    # -----------------------------------------------------------------------
                    # Switch to the next face of the cube
                    break

        # -----------------------------------------------------------------------
        # Sorting points along line of sight (there are at least two)
        # Arrange them in ascending and descending order.
        for alti_layer in range(nbi):
            for next_alti_layer in range(alti_layer + 1, nbi):
                if alti_layer_i[next_alti_layer] < alti_layer_i[alti_layer]:
                    dtmp = coord_col_i[alti_layer]
                    coord_col_i[alti_layer] = coord_col_i[next_alti_layer]
                    coord_col_i[next_alti_layer] = dtmp
                    dtmp = coord_row_i[alti_layer]
                    coord_row_i[alti_layer] = coord_row_i[next_alti_layer]
                    coord_row_i[next_alti_layer] = dtmp
                    dtmp = coord_alt_i[alti_layer]
                    coord_alt_i[alti_layer] = coord_alt_i[next_alti_layer]
                    coord_alt_i[next_alti_layer] = dtmp
                    dtmp = alti_layer_i[alti_layer]
                    alti_layer_i[alti_layer] = alti_layer_i[next_alti_layer]
                    alti_layer_i[next_alti_layer] = dtmp

        # -----------------------------------------------------------------------
        # Filtering points not located on the cube
        alti_layer = 0
        while alti_layer < nbi:
            # test inside the cube
            test_on_cube = (
                (coord_col_i[alti_layer] >= self.plane_coef_d[0])
                and (coord_col_i[alti_layer] <= self.plane_coef_d[1])
                and (coord_row_i[alti_layer] >= self.plane_coef_d[2])
                and (coord_row_i[alti_layer] <= self.plane_coef_d[3])
                and (coord_alt_i[alti_layer] >= self.plane_coef_d[4])
                and (coord_alt_i[alti_layer] <= self.plane_coef_d[5])
            )
            if not test_on_cube:
                # We translate all the following points (we overwrite this invalid point)
                for next_alti_layer in range(alti_layer + 1, nbi):
                    coord_col_i[next_alti_layer - 1] = coord_col_i[next_alti_layer]
                    coord_row_i[next_alti_layer - 1] = coord_row_i[next_alti_layer]
                    coord_alt_i[next_alti_layer - 1] = coord_alt_i[next_alti_layer]
                    alti_layer_i[next_alti_layer - 1] = alti_layer_i[next_alti_layer]
                nbi -= 1
            else:
                alti_layer += 1
        # -----------------------------------------------------------------------
        # No solution if 0 or 1 single point is found (we have tangent to the cube)
        if nbi < 2:
            b_trouve = False
            return True, b_trouve, point_b, h_intersect
        # -----------------------------------------------------------------------
        # There are only 2 points left so we cross the cube
        # LAIG-FA-MAJA-2168-CNES: no more filtering on identical points. There may be a number of points > 2
        # Init the current point
        # DTM coordinates
        point_dtm = np.zeros(3)
        point_dtm[0] = coord_col_i[0]
        point_dtm[1] = coord_row_i[0]
        point_dtm[2] = coord_alt_i[0]
        # point_dtm is the first intersection with the cube (line, column)
        # -----------------------------------------------------------------------
        # h is gld 3D
        h_intersect = alti_layer_i[0]
        # h_intersect is the h interpolation index (not integer)
        # -----------------------------------------------------------------------
        # Terrain coordinates
        point_b = self.index_to_ter(point_dtm)
        # point_b is the terrain point (lon,lat)
        # -----------------------------------------------------------------------
        # End, return
        b_trouve = True
        return (True, b_trouve, point_b, h_intersect)

    # gitlab issue #56
    # pylint: disable=too-many-locals
    # pylint: disable=too-many-function-args
    # pylint: disable=too-many-nested-blocks
    # pylint: disable=too-many-statements
    def intersection(self, los, point_b, h_intersect):  # noqa: C901
        """
        DTM intersection

        :param los:  line of sight
        :type los: numpy.array
        :param point_b:  position of intersection in DTM cube
        :type point_b: numpy.array
        :param h_intersect:  altitude in DTM cube
        :type h_intersect: float
        :return: intersection information (True,an intersection has been found ?, position of intersection)
        :rtype: tuple (bool, bool, numpy.array)
        """
        los_index = self.ters_to_indexs(los)
        point_b_dtm = self.ter_to_index(point_b)
        point_r = np.zeros(3)
        (npl, _) = los.shape
        alti = np.arange(npl, -1.0, -1.0)
        p_1 = point_b_dtm.copy()  # [p_1[0],p_1[1],p_1[2]]

        h_intersect_p1 = h_intersect

        n_row = self.dtm_image.nb_rows
        n_col = self.dtm_image.nb_columns

        # 1 - Init and preliminary tests
        #   1.1 - Test if the vertex is above the DTM
        #       - Compute DTM altitude ? vertex position
        alti_1 = self.interpolate(p_1[0], p_1[1])
        #       - Compute the altitude difference to the DTM
        d_alti_1 = p_1[2] - alti_1

        #       - Test if the new top point is above the DTM
        if d_alti_1 < 0:
            #       - The Point is below the DTM
            #          . means that the line of sight goes into the DTM by the side
            #          . then below, no solution.
            b_trouve = False
            return True, b_trouve, point_r

        #   1.2 - Init the rank of the first vertex of the line of sight
        i_0 = int(np.floor(h_intersect_p1))

        #   1.3 - Init the starting point (in p_2)
        p_2 = point_b_dtm.copy()
        h_intersect_p2 = h_intersect

        nb_planes = los_index.shape[0]
        # 2. - Loop on the grid planes
        while i_0 < (nb_planes - 1):
            # 2.1 - Init current vertex of line of sight
            col_0 = los_index[i_0][0]
            row_0 = los_index[i_0][1]
            z_0 = los_index[i_0][2]
            z_1 = los_index[i_0 + 1][2]

            # 2.2 - Init line of sight DTM
            los_dtm = los_index[i_0 + 1] - los_index[i_0]

            # 2.3 - Test if line of sight is vertical
            if los_dtm[0] == 0 and los_dtm[1] == 0:
                # 2.3.1 - LOS is  vertical:
                #    - Compute DTM altitude ? vertex position
                alti_1 = self.interpolate(col_0, row_0)

                #    Test if the next plane is above DTM
                if los_index[i_0 + 1][2] <= alti_1:
                    # Init exit point
                    p_1[0] = col_0
                    p_1[1] = row_0
                    p_1[2] = alti_1
                    b_trouve = True
                    point_r = self.index_to_ter(p_1)
                    return True, b_trouve, point_r
                # Positioning on next vertex
                i_0 += 1
            else:
                # 2.3.2 - LOS is not vertical :
                #         it will fly over the DTM
                #         . we can continue
                #         . remains to demonstrate that the LOS will crash on DTM...
                #
                # Init starting point
                # Init its LOS x-axis
                a_2 = h_intersect_p2 - i_0

                # Fixed an FA DG 10 bug, when  resetting, a_2 value can be 1
                # Which has the consequence of jumping to the next segment of the LOS
                # Thus, one can lose points on a slice of altitude
                # To be sure to scan the aiming segment, we reset
                # the starting point of the segment, i.e. a_2 = 0
                if a_2 >= 1.0:
                    a_2 = 0.0

                # Init first intersected mesh
                #  - Init mesh index
                col_c = int(np.floor(p_2[0]))
                row_c = int(np.floor(p_2[1]))

                # NB :    caution before to start:
                #        . we put ourselves on the right side of the mesh
                #        . in principle, you should not leave the DTM
                # We enter from the bottom, the DTM mesh is the previous one
                if (p_2[0] == col_c) and (los_dtm[0] < 0):
                    col_c -= 1

                # We enter from the left, the DTM mesh is the previous one
                if (p_2[1] == row_c) and (los_dtm[1] < 0):
                    row_c -= 1

                # LDD - We're already out of bounds, we stop
                if not ((a_2 < 1) and -1 < col_c < (n_row - 1) and -1 < row_c < (n_col - 1)):
                    b_trouve = False
                    return True, b_trouve, point_r

                # Iterative search loop of the intersected cell
                while (a_2 < 1) and -1 < col_c < (n_row - 1) and -1 < row_c < (n_col - 1):
                    # - Min and max altitudes of the mesh
                    h_i = self.alt_min_cell[col_c, row_c]
                    h_s = self.alt_max_cell[col_c, row_c]

                    # - Transfer: the low point becomes the high point
                    # a1 = a_2;
                    # p_1 becomes p_2
                    p_1 = p_2.copy()
                    h_intersect_p1 = h_intersect_p2

                    # 4.2 - Determination of a new low point
                    #      - LOS orientation test
                    if not los_dtm[0]:
                        # 4.2.1 - LOS is completely oriented  east-west
                        #   p_2[0] = p_1[0] ; // useless, is already init
                        #       - LOS orientation test
                        if los_dtm[1] < 0:
                            # 4.2.1.1 - LOS goes due west
                            p_2[1] = row_c
                            row_c -= 1
                        else:
                            # 4.2.1.2 - LOS goes due east
                            row_c += 1
                            p_2[1] = row_c

                        a_2 = (p_2[1] - row_0) / los_dtm[1]
                        p_2[2] = z_0 + a_2 * los_dtm[2]

                    elif not los_dtm[1]:
                        # 4.2.2 - LOS is oriented north-south
                        #  p_2[1] = p_1[1] ;
                        #       - LOS orientation test
                        if los_dtm[0] < 0:
                            # 4.2.2.1 - LOS goes due north
                            p_2[0] = col_c
                            col_c -= 1
                        else:
                            # 4.2.2.2 - LOS goes due south
                            col_c += 1
                            p_2[0] = col_c

                        a_2 = (p_2[0] - col_0) / los_dtm[0]
                        p_2[2] = z_0 + a_2 * los_dtm[2]
                    else:
                        # 4.2.3 - Any other LOS here
                        #            - Determination of the exit side
                        if (los_dtm[0] < 0) and (los_dtm[0] <= los_dtm[1]) and (los_dtm[0] <= -los_dtm[1]):
                            # 4.2.3.1 - LOS is mainly oriented north
                            #             - Intersect with north side
                            a_2 = (col_c - col_0) / los_dtm[0]
                            p_2[1] = row_0 + a_2 * los_dtm[1]

                            if (p_2[1] > row_c) and (p_2[1] < (row_c + 1)):
                                # LOS goes out by the north
                                p_2[0] = col_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[1] < row_c:
                                # LOS goes out by the west
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[1] > (row_c + 1):
                                # LOS goes out by the east
                                row_c += 1
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[1] == row_c:
                                # LOS goes out by the north-west corner
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1
                                row_c -= 1

                            elif p_2[1] == (row_c + 1):
                                # LOS goes out by the north-east corner
                                p_2[0] = col_c
                                row_c += 1
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                        elif (los_dtm[1] > 0) and (los_dtm[1] >= los_dtm[0]) and (los_dtm[1] >= -los_dtm[0]):
                            # 4.2.3.2 - LOS is mainly oriented east
                            #         - Intersect with east side
                            a_2 = (row_c + 1 - row_0) / los_dtm[1]
                            p_2[0] = col_0 + a_2 * los_dtm[0]

                            if (p_2[0] > col_c) and (p_2[0] < (col_c + 1)):
                                #  LOS goes out by the east
                                row_c += 1
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[0] < col_c:
                                # LOS goes out by the north
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[0] > (col_c + 1):
                                # LOS goes out by the south
                                col_c += 1
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[0] == col_c:
                                # LOS goes out by the north-east corner
                                row_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[0] == (col_c + 1):
                                # LOS goes out by the south-east corner
                                col_c += 1
                                row_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                        elif (los_dtm[0] > 0) and (los_dtm[0] >= los_dtm[1]) and (los_dtm[0] >= -los_dtm[1]):
                            # 4.2.3.3 - LOS is mainly oriented south
                            #         - Intersect with south side
                            a_2 = (col_c + 1 - col_0) / los_dtm[0]
                            p_2[1] = row_0 + a_2 * los_dtm[1]

                            if (p_2[1] > row_c) and (p_2[1] < (row_c + 1)):
                                # LOS goes out by the south
                                col_c += 1
                                p_2[0] = col_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[1] < row_c:
                                # LOS goes out by the west
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[1] > row_c + 1:
                                # LOS goes out by the east
                                row_c += 1
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[1] == row_c:
                                # LOS goes out by the south-west corner
                                col_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[1] == row_c + 1:
                                # LOS goes out by the south-east corner
                                col_c += 1
                                row_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                        elif (los_dtm[1] < 0) and (los_dtm[1] <= los_dtm[0]) and (los_dtm[1] <= -los_dtm[0]):
                            #  4.2.3.4 - VLOS is mainly oriented west
                            #          - Intersect with west side
                            a_2 = (row_c - row_0) / los_dtm[1]
                            p_2[0] = col_0 + a_2 * los_dtm[0]

                            if (p_2[0] > col_c) and (p_2[0] < col_c + 1):
                                #  LOS goes out by the west
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[0] < col_c:
                                #  LOS goes out by the north
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[0] > (col_c + 1):
                                #  LOS goes out by the south
                                col_c += 1
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[0] == col_c:
                                #  LOS goes out by the north-west corner
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1
                                row_c -= 1

                            elif p_2[0] == (col_c + 1):
                                #  LOS goes out by the south-west corner
                                col_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                    # LDD - min and max bounds of the "layer" Checking
                    b_intersect = False

                    if p_2[2] > z_0:
                        # We've gone too high, and that's not good!!!
                        b_intersect = not (((p_1[2] > h_s) and (z_0 > h_s)) or ((p_1[2] < h_i) and (z_0 < h_i)))

                    elif p_2[2] < z_1:
                        # We went too low, and that's not good either!!! (even if it already makes more sense)
                        b_intersect = not (((p_1[2] > h_s) and (z_1 > h_s)) or ((p_1[2] < h_i) and (z_1 < h_i)))

                    else:
                        b_intersect = not (((p_1[2] > h_s) and (p_2[2] > h_s)) or ((p_1[2] < h_i) and (p_2[2] < h_i)))

                    # 5. LOS intersection test with the cube

                    if b_intersect:
                        # There is intersection between LOS and the cube
                        # 5.1 - DTM Altitudes
                        alti_1 = self.interpolate(p_1[0], p_1[1])
                        h_2 = self.interpolate(p_2[0], p_2[1])

                        # 5.2 - Altitude differences with DTM
                        d_alti_1 = p_1[2] - alti_1
                        d_2 = p_2[2] - h_2

                        # 5.3 - Intersection test with DTM
                        if d_alti_1 * d_2 <= 0:
                            # There is intersection between los and the DTM
                            # 5.3.1 - Compute of approximate solution
                            d_2 = 2 * self.tol_z  # Init d_2 > TOL_Z
                            col_a = p_2[0]
                            row_a = p_2[1]
                            z_a = h_2

                            while abs(d_2) > self.tol_z:
                                # 5.3.1.1 - Linear interpolation coefficient of H
                                c_h = (p_1[2] - alti_1) / ((h_2 - alti_1) - (p_2[2] - p_1[2]))

                                # 5.3.1.2 - position of the interpolated point
                                col_a = p_1[0] + c_h * (p_2[0] - p_1[0])
                                row_a = p_1[1] + c_h * (p_2[1] - p_1[1])
                                z_a = p_1[2] + c_h * (p_2[2] - p_1[2])

                                # 5.3.1.3 - Altitude of the interpolated point
                                z_v = self.interpolate(col_a, row_a)

                                # 5.3.1.4 - Altitude difference of the interpolated point
                                d_2 = z_v - z_a

                                # 5.3.1.5 - Update
                                if d_2 < 0:
                                    # Update of the top point
                                    p_1[0] = col_a
                                    p_1[1] = row_a
                                    p_1[2] = z_a
                                    alti_1 = z_v

                                else:
                                    # Update of the low point
                                    p_2[0] = col_a
                                    p_2[1] = row_a
                                    p_2[2] = z_a
                                    h_2 = z_v

                            # // End, return
                            p_1[0] = col_a
                            p_1[1] = row_a
                            p_1[2] = z_a

                            b_trouve = True
                            point_r = self.index_to_ter(p_1)
                            return True, b_trouve, point_r

                # End loop on meshes

                # Test if we are still in the DTM cube si on est toujours dans le cube DTM
                if a_2 >= 1:
                    logging.debug("Change of plane")
                    # Change of plane
                    i_0 += 1

                    # Loading into p_2 of the new vertex
                    p_2 = los_index[i_0].copy()
                    h_intersect_p2 = alti[i_0]

                else:
                    # LDD - We looped on the meshes, we found nothing and we did not reach the next plane
                    # It means we're getting out of the grip, no need to continue
                    b_trouve = False
                    return True, b_trouve, point_r
            # End of general case (LOS not vertical)

        # End loop on the vertices
        # End, return
        b_trouve = False
        return True, b_trouve, point_r
