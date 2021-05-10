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
This module contains the euclidium help function
"""

import logging

DATUM = {"Z-M": "ellipsoid", "H-M": "geoid"}
CRS = {"WGS84": 4326, "GRS80": 4269}


def identify_gdlib_code(gdlib_code, default_datum="ellipsoid"):
    """
    Convert gdlib code to EPSG and datum.
    :param gdlib_code: gdlib string
    :type gdlib_code: str
    :param default_datum: default datum if gdlib code can't be parsed
    :type default_datum: str
    :returns: epsg code and datum
    :rtype: int,str
    """
    gdlib_code = gdlib_code.replace("/", "")
    splitted_code = gdlib_code.split(":")
    gdlib_srs = splitted_code[0]
    gdlib_datum = splitted_code[-1]
    datum = DATUM.get(gdlib_datum)
    if datum is None:
        logging.warning("gdlib datum %s is not correct default datum is set to ellipsoid", gdlib_datum)
        datum = default_datum
    if gdlib_srs.isdigit():
        epsg = int(gdlib_srs)
    else:
        epsg = CRS.get(gdlib_srs)
        if epsg is None:
            logging.warning("gdlib srs %s is not correct defaut srs is set to 4326", gdlib_srs)
            epsg = 4326
    return epsg, datum
