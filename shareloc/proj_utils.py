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
import osgeo
from osgeo import osr


def coordinates_conversion(coords, epsg_in, epsg_out):
    """
    Convert coords from a SRS to another one.
    :param coords: coords to project
    :type coords: numpy array
    :param epsg_in: EPSG code of the input SRS
    :type epsg_in: int
    :param epsg_out: EPSG code of the output SRS
    :type epsg_out: int
    :returns: converted coordinates
    :rtype: numpy array
    """
    srs_in = osr.SpatialReference()
    srs_in.ImportFromEPSG(epsg_in)
    srs_out = osr.SpatialReference()
    srs_out.ImportFromEPSG(epsg_out)
    # GDAL 3.0 Coordinate transformation (backwards compatibility)
    # https://github.com/OSGeo/gdal/issues/1546
    if int(osgeo.__version__[0]) >= 3:
        srs_in.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        srs_out.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    conversion = osr.CoordinateTransformation(srs_in, srs_out)
    coords = conversion.TransformPoints(coords)
    coords = np.array(coords)
    return coords
