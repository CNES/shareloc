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
This module contains the projection functions for shareloc
"""

import numpy as np
from rasterio import crs, warp


def coordinates_conversion(coords, epsg_in, epsg_out):
    """
    Convert coords from a SRS to another one.
    :param coords: coords to project
    :type coords: numpy array of 2D coords  (shape  (2,) or (N,2) or 3D coords (shape  (3,) or (N,3))
    :param epsg_in: EPSG code of the input SRS
    :type epsg_in: int
    :param epsg_out: EPSG code of the output SRS
    :type epsg_out: int
    :returns: converted coordinates
    :rtype: numpy array of 2D coord (N,2) or 3D coords (N,3)
    """
    srs_in = crs.CRS.from_epsg(epsg_in)
    srs_out = crs.CRS.from_epsg(epsg_out)
    if (coords.size / 3 == 1 or coords.size / 2 == 1) and (coords.ndim == 1):
        coords = coords[np.newaxis, :]
    alti = None
    if coords.shape[1] == 3:
        alti = coords[:, 2]
    coords = np.array(warp.transform(srs_in, srs_out, coords[:, 0], coords[:, 1], alti))
    coords = coords.transpose()
    return coords
