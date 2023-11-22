#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2022 Centre National d'Etudes Spatiales (CNES).
# Copyright (c) 2023 CS GROUP - France, https://csgroup.eu
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
from os.path import basename
from typing import Dict
from xml.dom import minidom

# Third party imports
import numpy as np
import rasterio as rio
from numba import config, njit, prange

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.proj_utils import coordinates_conversion


import sys
sys.path.append(".")
import libs.pbrpc as bind


# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"



@GeoModel.register("RPC_optim")
class RPC_optim(bind.RPC,GeoModelTemplate):
    """
    RPC optimizes with cpp bindings class including direct and inverse localization instance methods
    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self,geomodel_path: str):

        bind.RPC.__init__(self)
        GeoModelTemplate.__init__(self,geomodel_path)

        self.type = "RPC_optim"

        # RPC parameters are load from geomodel_path to rpc params
        self.load()