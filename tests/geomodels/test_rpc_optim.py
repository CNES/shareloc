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


import os

import numpy as np

# Third party imports
import pytest

# Shareloc bindings
import rpc_c

# Shareloc imports
from shareloc.geomodels import GeoModel
from shareloc.geomodels.rpc import compute_rational_function_polynomial, polynomial_equation
from shareloc.geomodels.rpc_readers import rpc_reader

# Shareloc test imports
from ..helpers import data_path


# pylint: disable=duplicate-code
@pytest.mark.parametrize(
    "geom_path",
    [
        "rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom",
        "rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif",
        "rpc/PHRDIMAP_P1BP--2017030824934340CP.XML",
        "rpc/RPC_P1BP--2017092838284574CP.XML",
        "rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML",
    ],
)
def test_construtor(geom_path):
    """
    Test RpcOptim constructor
    """

    file_path = os.path.join(data_path(), geom_path)

    rpc_py = GeoModel(file_path, "RPC")

    rpc_params = rpc_reader(file_path, topleftconvention=True)
    norm_coeffs = [
        rpc_params["offset_x"],
        rpc_params["scale_x"],  # longitude
        rpc_params["offset_y"],
        rpc_params["scale_y"],  # latitude
        rpc_params["offset_alt"],
        rpc_params["scale_alt"],
        rpc_params["offset_col"],
        rpc_params["scale_col"],
        rpc_params["offset_row"],
        rpc_params["scale_row"],
    ]
    rpc_cpp = rpc_c.RPC(
        rpc_params["num_col"], rpc_params["den_col"], rpc_params["num_row"], rpc_params["den_row"], norm_coeffs
    )

    assert rpc_py.offset_x == rpc_cpp.get_offset_lon()
    assert rpc_py.scale_x == rpc_cpp.get_scale_lon()
    assert rpc_py.offset_y == rpc_cpp.get_offset_lat()
    assert rpc_py.scale_y == rpc_cpp.get_scale_lat()
    assert rpc_py.offset_alt == rpc_cpp.get_offset_alt()
    assert rpc_py.scale_alt == rpc_cpp.get_scale_alt()
    assert rpc_py.offset_col == rpc_cpp.get_offset_col()
    assert rpc_py.scale_col == rpc_cpp.get_scale_col()
    assert rpc_py.offset_row == rpc_cpp.get_offset_row()
    assert rpc_py.scale_row == rpc_cpp.get_scale_row()

    np.testing.assert_array_equal(rpc_py.num_col, rpc_cpp.get_num_col())
    np.testing.assert_array_equal(rpc_py.den_col, rpc_cpp.get_den_col())
    np.testing.assert_array_equal(rpc_py.num_row, rpc_cpp.get_num_row())
    np.testing.assert_array_equal(rpc_py.den_row, rpc_cpp.get_den_row())


def test_method_rpc_cpp():
    """
    Test call to method parent class rpc in cpp
    """

    file_path = os.path.join(data_path(), "rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom")
    rpc = GeoModel(file_path, "RpcOptim")

    vector_double = [1.0, 1.0, 1.0]
    double = 1.0
    integer = 1
    string = "string"

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
    array20 = [1.0 for i in range(20)]
    double = 1.0

    rpc_c.polynomial_equation(double, double, double, array20)

    rpc_c.compute_rational_function_polynomial(
        vector_double,
        vector_double,
        vector_double,
        array20,
        array20,
        array20,
        array20,
        double,
        double,
        double,
        double,
    )

    rpc_c.derivative_polynomial_latitude(double, double, double, vector_double)

    rpc_c.derivative_polynomial_longitude(double, double, double, vector_double)

    rpc_c.compute_loc_inverse_derivates_optimized(
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


@pytest.mark.parametrize(
    "geom_path",
    [
        "rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom",
        "rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif",
        "rpc/PHRDIMAP_P1BP--2017030824934340CP.XML",
        "rpc/RPC_P1BP--2017092838284574CP.XML",
        "rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML",
    ],
)
def test_polynomial_equation(geom_path):
    """
    test polynomial_equation function
    """

    file_path = os.path.join(data_path(), geom_path)
    rpc_params = rpc_reader(file_path, topleftconvention=True)
    norm_coeffs = [
        rpc_params["offset_x"],
        rpc_params["scale_x"],  # longitude
        rpc_params["offset_y"],
        rpc_params["scale_y"],  # latitude
        rpc_params["offset_alt"],
        rpc_params["scale_alt"],
        rpc_params["offset_col"],
        rpc_params["scale_col"],
        rpc_params["offset_row"],
        rpc_params["scale_row"],
    ]
    rpc_cpp = rpc_c.RPC(
        rpc_params["num_col"], rpc_params["den_col"], rpc_params["num_row"], rpc_params["den_row"], norm_coeffs
    )

    # arbitrary values (extract from a rpc.py test)
    xnorm = -0.95821893  # lon_norm
    ynorm = 0.97766941  # lat_norm
    znorm = -5.29411765  # alt_norm

    res_c_den_col = rpc_c.polynomial_equation(xnorm, ynorm, znorm, rpc_cpp.get_den_col())
    res_c_den_row = rpc_c.polynomial_equation(xnorm, ynorm, znorm, rpc_cpp.get_den_row())
    res_c_num_col = rpc_c.polynomial_equation(xnorm, ynorm, znorm, rpc_cpp.get_num_col())
    res_c_num_row = rpc_c.polynomial_equation(xnorm, ynorm, znorm, rpc_cpp.get_num_row())

    # rpc.py polynomial_equation
    res_py_den_col = polynomial_equation(xnorm, ynorm, znorm, np.array(rpc_cpp.get_den_col(), dtype=np.float64))
    res_py_den_row = polynomial_equation(xnorm, ynorm, znorm, np.array(rpc_cpp.get_den_row(), dtype=np.float64))
    res_py_num_col = polynomial_equation(xnorm, ynorm, znorm, np.array(rpc_cpp.get_num_col(), dtype=np.float64))
    res_py_num_row = polynomial_equation(xnorm, ynorm, znorm, np.array(rpc_cpp.get_num_row(), dtype=np.float64))

    assert res_c_den_col == pytest.approx(res_py_den_col, abs=1e-15)
    assert res_c_den_row == pytest.approx(res_py_den_row, abs=1e-15)
    assert res_c_num_col == pytest.approx(res_py_num_col, abs=1e-15)
    assert res_c_num_row == pytest.approx(res_py_num_row, abs=1e-15)


@pytest.mark.parametrize(
    "geom_path",
    [
        "rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom",
        "rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif",
        "rpc/PHRDIMAP_P1BP--2017030824934340CP.XML",
        "rpc/RPC_P1BP--2017092838284574CP.XML",
        "rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML",
    ],
)
def test_compute_rational_function_polynomial(geom_path):
    """
    test compute_rational_function_polynomial function
    """

    file_path = os.path.join(data_path(), geom_path)
    rpc_params = rpc_reader(file_path, topleftconvention=True)
    norm_coeffs = [
        rpc_params["offset_x"],
        rpc_params["scale_x"],  # longitude
        rpc_params["offset_y"],
        rpc_params["scale_y"],  # latitude
        rpc_params["offset_alt"],
        rpc_params["scale_alt"],
        rpc_params["offset_col"],
        rpc_params["scale_col"],
        rpc_params["offset_row"],
        rpc_params["scale_row"],
    ]
    rpc_cpp = rpc_c.RPC(
        rpc_params["num_col"], rpc_params["den_col"], rpc_params["num_row"], rpc_params["den_row"], norm_coeffs
    )

    # arbitrary values (extract from a rpc.py test) + last different one
    xnorm = [-0.95821893, 0.0, -0.96038983, -0.95821891, -0.95821893, 0.354]  # lon_norm
    ynorm = [0.97766941, 0.0, 0.98172072, 0.97766945, 0.97766941, -0.654]  # lat_norm
    znorm = [-5.29411765, -5.29411765, -5.29411765, -5.29411765, -5.29411765, 5.1]  # alt_norm

    res_cpp = rpc_c.compute_rational_function_polynomial(
        xnorm,
        ynorm,
        znorm,
        rpc_cpp.get_num_col(),
        rpc_cpp.get_den_col(),
        rpc_cpp.get_num_row(),
        rpc_cpp.get_den_row(),
        rpc_cpp.get_scale_col(),
        rpc_cpp.get_offset_col(),
        rpc_cpp.get_scale_row(),
        rpc_cpp.get_offset_row(),
    )

    res_py = compute_rational_function_polynomial(
        np.array(xnorm, dtype=np.float64),
        np.array(ynorm, dtype=np.float64),
        np.array(znorm, dtype=np.float64),
        np.array(rpc_cpp.get_num_col(), dtype=np.float64),
        np.array(rpc_cpp.get_den_col(), dtype=np.float64),
        np.array(rpc_cpp.get_num_row(), dtype=np.float64),
        np.array(rpc_cpp.get_den_row(), dtype=np.float64),
        rpc_cpp.get_scale_col(),
        rpc_cpp.get_offset_col(),
        rpc_cpp.get_scale_row(),
        rpc_cpp.get_offset_row(),
    )

    np.testing.assert_array_almost_equal(np.array(res_cpp[0]), res_py[0], 11)
    np.testing.assert_array_almost_equal(np.array(res_cpp[0]), res_py[0], 11)
