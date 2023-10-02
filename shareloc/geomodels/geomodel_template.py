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
from typing import Dict

# Third party imports


# pylint: disable=too-few-public-methods
class GeoModelTemplate(metaclass=ABCMeta):
    """
    Class for general specification of a geometric model
    declined in rpc.py and grid.py
    """

    @abstractmethod
    def __init__(self, geomodel_file: str, geomodel_params: Dict = None):
        """
        Return the geomodel object associated with the geomodel_type
        given in the configuration

        :param geomodel_file: path of the geomodel file to instanciate.
        :type geomodel_file: string
        :param geomodel_params: Geomodels parameters as dict (unused, just to fit RPC) TODO: evolve with new API
        :type geomodel_params: Dict
        """
        # geomodel filename path
        self.filename = geomodel_file

        # geomodel_params if exists (to evolve)
        self.geomodels_params = geomodel_params

        #
