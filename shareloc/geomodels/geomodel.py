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
This module contains the GeoModel class factory.
This main API GeoModel class generates an object linked
with geomodel configuration "geomodel_name"
(from registered GeoModelTemplate class)

** Experimental and not used as API for now**
** Will be used with an automatic switcher between Grid file format and
RPC file format obj = GeoModel("file_geom") only without geomodel_type
"""

# Standard imports
import logging
from typing import Any, Dict


class GeoModel:
    """
    GeoModel factory:
    A class designed for registered all available geomodels
    and instantiate them when needed.
    """

    # Dict (geomodel_name: str, class: object) containing registered geomodels
    available_geomodels: Dict[str, Any] = {}

    def __new__(cls, geomodel_path: str, geomodel_type: str = "RPC"):
        """
        Return a GeoModelTemplate child instance
        associated with the "geomodel_name" given in the configuration
        through create_geomodel local method for clarity.

        TODO: optional geomodel_type would profit to have an automatic switcher between geomodels (RPC, grids, ...)

        :param geomodel_path: Path of geomodel file
        :type geomodel_path: string
        :param geomodel_type: Type of the geomodel, default "RPC", used as key for specific geomodel subclass instance
        :type geomodel_type: string
        """
        return cls.create_geomodel(geomodel_path, geomodel_type)

    @classmethod
    def create_geomodel(cls, geomodel_path: str, geomodel_type: str = "RPC"):
        """
        Factory command to create the geomodel from geomodel_type
        Return a GeoModelTemplate child instance
        associated with the "geomodel_type"

        :param geomodel_path: Path of geomodel file
        :type geomodel_path: string
        :param geomodel_type: Type of the geomodel, default "RPC", used as key for specific geomodel subclass instance
        :type geomodel_type: string
        """
        # Create geomodel object with geomodel_path parameter from geomodel_type if exists
        try:
            geomodel_class = cls.available_geomodels[geomodel_type]
            geomodel = geomodel_class.load(geomodel_path)
            logging.debug("GeoModel type name: %s", geomodel_type)
        except KeyError:
            logging.error("Geomodel type named %s is not supported", geomodel_type)
            raise

        return geomodel

    @classmethod
    def print_avalaible_geomodels(cls):
        """
        Print all registered applications
        """
        for geomodel_type in cls.available_geomodels:
            print(geomodel_type)

    @classmethod
    def register(cls, geomodel_type: str):
        """
        Allows to register the GeoModelTemplate subclass in the factory
        with its geomodel geomodel_type through decorator

        :param geomodel_type: the subclass name to be registered
        :type geomodel_type: string
        """

        def decorator(geomodel_subclass):
            """
            Register the geomodel subclass in the available methods

            :param geomodel_subclass: the subclass to be registered
            :type geomodel_subclass: object
            """
            cls.available_geomodels[geomodel_type] = geomodel_subclass
            return geomodel_subclass

        return decorator
