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
Module to test RpcOptim class
"""

# Third party imports
import numpy as np
import pytest

import rpc_c

# Shareloc imports
from shareloc.geomodels import GeoModel


@pytest.mark.parametrize(
    "geom_path",
    [
        "tests/data/rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom",
        "tests/data/rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif",
        "tests/data/rpc/PHRDIMAP_P1BP--2017030824934340CP.XML",
        "tests/data/rpc/RPC_P1BP--2017092838284574CP.XML",
        "tests/data/rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML",
    ],
)
def test_load_rpc_params(geom_path):
    """
    test loading of rpc_params
    """
    geom = GeoModel(geom_path, "RpcOptim").__dict__
    geom_ref = GeoModel(geom_path).__dict__

    del geom["type"]
    del geom_ref["type"]

    for key, value in geom.items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(value, geom_ref[key])
        else:
            assert value == geom_ref[key]


def test_method_rpc_cpp():
    """
    Test call to method parent class rpc in cpp
    """

    rpc = GeoModel("tests/data/rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom", "RpcOptim")

    vector_double = [1.0, 1.0, 1.0]
    double = 1.0
    integer = 1
    string = "peuimporte"

    rpc.direct_loc_h(vector_double, vector_double, double, False)
    rpc.direct_loc_grid_h(integer, integer, integer, integer, integer, integer, double)
    rpc.direct_loc_dtm(double, double, string)
    rpc.inverse_loc(vector_double, vector_double, double)
    rpc.filter_coordinates(vector_double, vector_double, False, string)
    rpc.compute_loc_inverse_derivates(vector_double, vector_double, vector_double)
    rpc.direct_loc_inverse_iterative(vector_double, vector_double, double, integer, False)
    rpc.get_alt_min_max()
    rpc.los_extrema(double, double, double, double, False, integer)


def test_function_rpc_cpp():
    """
    Test call to function written in rpc.cpp
    """

    vector_double = [1.0, 1.0, 1.0]
    double = 1.0

    rpc_c.polynomial_equation(double, double, double, vector_double)

    rpc_c.compute_rational_function_polynomial(
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        double,
        double,
        double,
        double,
    )

    rpc_c.derivative_polynomial_latitude(double, double, double, vector_double)

    rpc_c.derivative_polynomial_longitude(double, double, double, vector_double)

    rpc_c.compute_loc_inverse_derivates_numba(
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        vector_double,
        double,
        double,
        double,
        double,
    )
