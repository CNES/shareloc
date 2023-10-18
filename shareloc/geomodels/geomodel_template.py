#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of shareloc
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
This module contains the coregistration class template.
It contains the structure for all coregistration methods in subclasses and
generic coregistration code to avoid duplication.
"""

# Standard imports
from abc import ABCMeta, abstractmethod

# Global variable for optimization mode (functions in C)
# SHARELOC_OPTIM_GEOMODEL = False

# TODO: Override functions depending on optimization or not

# if(SHARELOC_OPTIM_GEOMODEL == True):
#     GeoModelTemplate.direct_loc_dtm = GeoModelTemplate.direct_loc_dtm_optim

# pylint: disable=too-few-public-methods


class GeoModelTemplate(metaclass=ABCMeta):
    """
    Class for general specification of a geometric model
    declined in rpc.py and grid.py
    """

    @abstractmethod
    def __init__(self, geomodel_path: str):
        """
        Return the geomodel object associated with the geomodel_type
        given in the configuration

        :param geomodel_path: path of the geomodel file to instanciate.
        :type geomodel_path: string
        """
        # geomodel filename path
        self.geomodel_path: str = geomodel_path

        # geomodel type. Set by the subclass
        self.type: str
