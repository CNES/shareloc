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
This module contains the optimized (with cpp bindings) RPC class corresponding to the RPC models.
RPC models covered are : DIMAP V1, DIMAP V2, DIMAP V3, ossim (geom file), geotiff.
"""


import numpy as np

# Third party imports
from numba import config

import rpc_c as bind

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.geomodels.rpc_readers import rpc_reader

# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"


# pylint: disable=duplicate-code
@GeoModel.register("RpcOptim")
class RpcOptim(bind.RPC, GeoModelTemplate):
    """
    RPC optimized with cpp bindings class including direct and inverse localization instance methods
    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self, rpc_params):
        bind.RPC.__init__(self)
        GeoModelTemplate.__init__(self)

        self.offset_alt = None
        self.scale_alt = None
        self.offset_col = None
        self.scale_col = None
        self.offset_row = None
        self.scale_row = None
        self.offset_x = None
        self.scale_x = None
        self.offset_y = None
        self.scale_y = None

        self.datum = None
        for key, value in rpc_params.items():
            setattr(self, key, value)

        self.type = "RpcOptim"
        if self.epsg is None:
            self.epsg = 4326
        if self.datum is None:
            self.datum = "ellipsoid"

        self.lim_extrapol = 1.0001

        # Each monome: c[0]*X**c[1]*Y**c[2]*Z**c[3]
        monomes_order = [
            [1, 0, 0, 0],
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [1, 0, 0, 1],
            [1, 1, 1, 0],
            [1, 1, 0, 1],
            [1, 0, 1, 1],
            [1, 2, 0, 0],
            [1, 0, 2, 0],
            [1, 0, 0, 2],
            [1, 1, 1, 1],
            [1, 3, 0, 0],
            [1, 1, 2, 0],
            [1, 1, 0, 2],
            [1, 2, 1, 0],
            [1, 0, 3, 0],
            [1, 0, 1, 2],
            [1, 2, 0, 1],
            [1, 0, 2, 1],
            [1, 0, 0, 3],
        ]

        self.monomes = np.array(monomes_order)

        # monomial coefficients of 1st variable derivative
        self.monomes_deriv_1 = np.array(
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [1, 0, 1, 0],
                [1, 0, 0, 1],
                [0, 0, 1, 1],
                [2, 1, 0, 0],
                [0, 0, 2, 0],
                [0, 0, 0, 2],
                [1, 0, 1, 1],
                [3, 2, 0, 0],
                [1, 0, 2, 0],
                [1, 0, 0, 2],
                [2, 1, 1, 0],
                [0, 0, 3, 0],
                [0, 0, 1, 2],
                [2, 1, 0, 1],
                [0, 0, 2, 1],
                [0, 0, 0, 3],
            ]
        )

        # monomial coefficients of 2nd variable derivative
        self.monomes_deriv_2 = np.array(
            [
                [0, 0, 0, 0],
                [0, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 0, 1],
                [1, 1, 0, 0],
                [0, 1, 0, 1],
                [1, 0, 0, 1],
                [0, 2, 0, 0],
                [2, 0, 1, 0],
                [0, 0, 0, 2],
                [1, 1, 0, 1],
                [0, 3, 0, 0],
                [2, 1, 1, 0],
                [0, 1, 0, 2],
                [1, 2, 0, 0],
                [3, 0, 2, 0],
                [1, 0, 0, 2],
                [0, 2, 0, 1],
                [2, 0, 1, 1],
                [0, 0, 0, 3],
            ]
        )

        self.inverse_coefficient = False
        self.direct_coefficient = False

        # pylint: disable=access-member-before-definition
        if self.num_col:
            self.inverse_coefficient = True
            self.num_col = np.array(self.num_col)
            self.den_col = np.array(self.den_col)
            self.num_row = np.array(self.num_row)
            self.den_row = np.array(self.den_row)

        # pylint: disable=access-member-before-definition
        if self.num_x:
            self.direct_coefficient = True
            self.num_x = np.array(self.num_x)
            self.den_x = np.array(self.den_x)
            self.num_y = np.array(self.num_y)
            self.den_y = np.array(self.den_y)

        self.alt_minmax = [self.offset_alt - self.scale_alt, self.offset_alt + self.scale_alt]
        self.col0 = self.offset_col - self.scale_col
        self.colmax = self.offset_col + self.scale_col
        self.row0 = self.offset_row - self.scale_row
        self.rowmax = self.offset_row + self.scale_row

    @classmethod
    def load(cls, geomodel_path):
        """
        Load from any RPC (auto identify driver)
        from filename (dimap, ossim kwl, geotiff)

        TODO: topleftconvention always to True, set a standard and remove the option

        topleftconvention boolean: [0,0] position
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        # Set topleftconvention (keeping historic option): to clean
        cls.geomodel_path = geomodel_path
        return cls(rpc_reader(geomodel_path, topleftconvention=True))
