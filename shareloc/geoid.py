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
This module contains the Geoid class to retrieve geoid height above ellispoid.
"""

import numpy as np
from scipy import interpolate
from shareloc.image import Image


def interpolate_geoid_height(geoid_filename, positions, interpolation_method="linear"):
    """
    terrain to index conversion
    :param geoid_filename: geoid_filename
    :type geoid_filename: str
    :param positions: geodetic coordinates
    :type positions: 2D numpy array : (number of points, [long coord, lat coord])
    :parama interpolation_method default is 'linear' (interpn interpolation method)
    :type str
    :return geoid height
    :rtype 1 numpy array (nuber of points)
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
        raise RuntimeError("Geoid can" "t handle latitudes greater than 90 deg.")
    indexes_geoid = geoid_image.transform_physical_point_to_index(positions[:, 1], positions[:, 0])
    return interpolate.interpn(points, geoid_image.data[:, :], indexes_geoid, method=interpolation_method)
