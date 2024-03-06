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

"""
Module to test RPCoptim class
"""


import os

import numpy as np

# Third party imports
import pytest

# Shareloc bindings
import bindings_cpp
from shareloc.dtm_reader import dtm_reader
from shareloc.geofunctions.dtm_intersection import DTMIntersection

# Shareloc imports
from shareloc.geomodels import GeoModel
from shareloc.geomodels.rpc import (
    compute_rational_function_polynomial,
    derivative_polynomial_latitude,
    derivative_polynomial_longitude,
    polynomial_equation,
)

# Shareloc test imports
from ..helpers import bindings_cpp_constructor, data_path


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
    Test RPCoptim constructor
    """

    file_path = os.path.join(data_path(), geom_path)

    rpc_py = GeoModel(file_path, "RPC")
    rpc_cpp = bindings_cpp_constructor(file_path)

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
    rpc_cpp = bindings_cpp_constructor(file_path)

    # arbitrary values (extract from a rpc.py test)
    xnorm = -0.95821893  # lon_norm
    ynorm = 0.97766941  # lat_norm
    znorm = -5.29411765  # alt_norm

    coord_array = [
        1.00000000000000000,
        -0.95821893000000002,
        0.97766940999999996,
        -5.29411764999999956,
        -0.93682133594393124,
        5.07292374987711447,
        -5.17589687934608556,
        0.91818351781034491,
        0.95583747524974805,
        28.02768169204151860,
        4.95964236951734527,
        -0.87982082797986472,
        -0.91590156278771506,
        -26.85665516132861441,
        0.89767993812936431,
        0.93449306048331071,
        27.40180702352602893,
        -4.86097156757883564,
        -5.06031604825112868,
        -148.38184433441884380,
    ]

    res_c_den_col = bindings_cpp.polynomial_equation(coord_array, rpc_cpp.get_den_col())
    res_c_den_row = bindings_cpp.polynomial_equation(coord_array, rpc_cpp.get_den_row())
    res_c_num_col = bindings_cpp.polynomial_equation(coord_array, rpc_cpp.get_num_col())
    res_c_num_row = bindings_cpp.polynomial_equation(coord_array, rpc_cpp.get_num_row())

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
    rpc_cpp = bindings_cpp_constructor(file_path)

    # arbitrary values (extract from a rpc.py test) + last different one
    xnorm = [-0.95821893, 0.0, -0.96038983, -0.95821891, -0.95821893, 0.354, np.nan, 15, 15]  # lon_norm len+1
    ynorm = [0.97766941, 0.0, 0.98172072, 0.97766945, 0.97766941, -0.654, 15, np.nan]  # lat_norm
    znorm = [-5.29411765, -5.29411765, -5.29411765, -5.29411765, -5.29411765, 5.1, 15, 15]  # alt_norm

    res_cpp = rpc_cpp.compute_rational_function_polynomial(
        xnorm,
        ynorm,
        znorm,
        rpc_cpp.get_num_col(),
        rpc_cpp.get_den_col(),
        rpc_cpp.get_num_row(),
        rpc_cpp.get_den_row(),
        False,
        "inverse",
        rpc_cpp.get_scale_lon(),
        rpc_cpp.get_offset_lon(),
        rpc_cpp.get_scale_lat(),
        rpc_cpp.get_offset_lat(),
        rpc_cpp.get_scale_alt(),
        rpc_cpp.get_offset_alt(),
        rpc_cpp.get_scale_col(),
        rpc_cpp.get_offset_col(),
        rpc_cpp.get_scale_row(),
        rpc_cpp.get_offset_row(),
    )
    # compute_rational_function_polynomial_unitary
    res_cpp_row = np.empty((len(znorm)))
    res_cpp_col = np.empty((len(znorm)))

    for i, znorm_i in enumerate(znorm):
        col_i, row_i, _ = rpc_cpp.compute_rational_function_polynomial_unitary(
            xnorm[i],
            ynorm[i],
            znorm_i,
            rpc_cpp.get_num_col(),
            rpc_cpp.get_den_col(),
            rpc_cpp.get_num_row(),
            rpc_cpp.get_den_row(),
            False,
            "inverse",
            rpc_cpp.get_scale_lon(),
            rpc_cpp.get_offset_lon(),
            rpc_cpp.get_scale_lat(),
            rpc_cpp.get_offset_lat(),
            rpc_cpp.get_scale_alt(),
            rpc_cpp.get_offset_alt(),
            rpc_cpp.get_scale_col(),
            rpc_cpp.get_offset_col(),
            rpc_cpp.get_scale_row(),
            rpc_cpp.get_offset_row(),
        )
        res_cpp_row[i] = row_i
        res_cpp_col[i] = col_i

    # Ref python
    xnorm = (np.array(xnorm[:-3]) - rpc_cpp.get_offset_lon()) / rpc_cpp.get_scale_lon()
    ynorm = (np.array(ynorm[:-2]) - rpc_cpp.get_offset_lat()) / rpc_cpp.get_scale_lat()
    znorm = (np.array(znorm[:-2]) - rpc_cpp.get_offset_alt()) / rpc_cpp.get_scale_alt()

    res_py = compute_rational_function_polynomial(
        xnorm,
        ynorm,
        znorm,
        np.array(rpc_cpp.get_num_col(), dtype=np.float64),
        np.array(rpc_cpp.get_den_col(), dtype=np.float64),
        np.array(rpc_cpp.get_num_row(), dtype=np.float64),
        np.array(rpc_cpp.get_den_row(), dtype=np.float64),
        rpc_cpp.get_scale_col(),
        rpc_cpp.get_offset_col(),
        rpc_cpp.get_scale_row(),
        rpc_cpp.get_offset_row(),
    )

    res_py_0 = np.append(res_py[0], [np.nan, np.nan])
    res_py_1 = np.append(res_py[1], [np.nan, np.nan])

    np.testing.assert_allclose(np.array(res_cpp[1]), res_py_0, 0, 3e-9)
    np.testing.assert_allclose(np.array(res_cpp[0]), res_py_1, 0, 2e-9)

    np.testing.assert_allclose(res_cpp_row, res_py_0, 0, 3e-9)
    np.testing.assert_allclose(res_cpp_col, res_py_1, 0, 2e-9)


@pytest.mark.parametrize(
    "id_scene,lon,lat,alt, col_vt,row_vt",
    [
        (
            "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.geom",
            [7.048662660737769592],
            [43.72774839443545858],
            [0.0],
            100.5,
            200.5,
        ),
        (
            "RPC_PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.XML",
            [7.048662660737769592],
            [43.72774839443545858],
            [0.0],
            100.5,
            200.5,
        ),
        (
            "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif",
            [7.048662660737769592],
            [43.72774839443545858],
            [0.0],
            100.5,
            200.5,
        ),
    ],
)
def test_inverse_loc_from_any_input(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization from any file and differents configuration of inputs
    """
    data_folder = data_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)

    rpc_optim = GeoModel(rpc_file, "RPCoptim")

    (row, col, alt_res) = rpc_optim.inverse_loc(lon, lat, alt)

    assert col[0] == pytest.approx(col_vt, abs=1e-2)
    assert row[0] == pytest.approx(row_vt, abs=1e-2)
    assert alt_res[0] == alt[0]

    # Check arg of differents sizes
    lon_ext = lon * 3
    lat_ext = lat * 2
    (row, col, alt_res) = rpc_optim.inverse_loc(lon_ext, lat_ext, alt)

    assert len(row) == 2
    assert len(alt_res) == 2
    assert len(alt_res) == 2
    assert col[0] == pytest.approx(col_vt, abs=1e-2)
    assert row[0] == pytest.approx(row_vt, abs=1e-2)
    assert alt_res[0] == alt[0]

    lon_ext = lon * 2
    lat_ext = lat * 3
    (row, col, alt_res) = rpc_optim.inverse_loc(lon_ext, lat_ext, alt)

    assert len(row) == 2
    assert len(alt_res) == 2
    assert len(alt_res) == 2
    assert col[0] == pytest.approx(col_vt, abs=1e-2)
    assert row[0] == pytest.approx(row_vt, abs=1e-2)
    assert alt_res[0] == alt[0]


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
def test_inverse_loc(geom_path):
    """
    test inverse localization accuracy
    """

    rpc_path = os.path.join(data_path(), geom_path)

    rpc_optim = GeoModel(rpc_path, "RPCoptim")
    rpc_py = GeoModel(rpc_path, "RPC")

    # INPUTS
    nrb_point = 1e6
    first_lon = 7.0477886581984
    first_lat = 43.62208491280199
    last_lon = 7.308411551163017
    last_lat = 43.73298365695963
    first_alt = -50
    last_alt = 1000

    lon_vect = np.linspace(first_lon, last_lon, int(nrb_point ** (1 / 3)) + 1)
    lat_vect = np.linspace(first_lat, last_lat, int(nrb_point ** (1 / 3)) + 1)
    alt_vect = np.linspace(first_alt, last_alt, int(nrb_point ** (1 / 3)) + 1)

    lon_vect, lat_vect, alt_vect = np.meshgrid(lon_vect, lat_vect, alt_vect)

    lon_vect = np.ndarray.flatten(lon_vect)
    lat_vect = np.ndarray.flatten(lat_vect)
    alt_vect = np.ndarray.flatten(alt_vect)

    lon_vect[0] = np.nan
    res_cpp = rpc_optim.inverse_loc(lon_vect, lat_vect, alt_vect)
    res_py = rpc_py.inverse_loc(lon_vect, lat_vect, alt_vect)

    np.testing.assert_allclose(res_cpp[0], res_py[0], 0, 2e-8)
    np.testing.assert_allclose(res_cpp[1], res_py[1], 0, 2e-8)

    # Inverse loc unitary

    file_path = os.path.join(data_path(), rpc_path)
    rpc_cpp = bindings_cpp_constructor(file_path)

    lon_out = np.empty((len(lon_vect)))
    lat_out = np.empty((len(lat_vect)))

    for i, lon_vect_i in enumerate(lon_vect):
        lon_i, lat_i, __ = rpc_cpp.inverse_loc(lon_vect_i, lat_vect[i], alt_vect[i])
        lon_out[i] = lon_i
        lat_out[i] = lat_i

    np.testing.assert_allclose(lon_out, res_py[0], 0, 2e-8)
    np.testing.assert_allclose(lat_out, res_py[1], 0, 2e-8)


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
def test_derivative_polynomial_latitude(geom_path):
    """
    test derivative_polynomial_latitude function
    """

    file_path = os.path.join(data_path(), geom_path)
    rpc_cpp = bindings_cpp_constructor(file_path)

    # arbitrary values (extract from a rpc.py test)
    xnorm = -0.95821893  # lon_norm
    ynorm = 0.97766941  # lat_norm
    znorm = -5.29411765  # alt_norm

    res_c_den_col = bindings_cpp.derivative_polynomial_latitude(xnorm, ynorm, znorm, rpc_cpp.get_den_col())
    res_c_den_row = bindings_cpp.derivative_polynomial_latitude(xnorm, ynorm, znorm, rpc_cpp.get_den_row())
    res_c_num_col = bindings_cpp.derivative_polynomial_latitude(xnorm, ynorm, znorm, rpc_cpp.get_num_col())
    res_c_num_row = bindings_cpp.derivative_polynomial_latitude(xnorm, ynorm, znorm, rpc_cpp.get_num_row())

    # rpc.py derivative_polynomial_latitude
    res_py_den_col = derivative_polynomial_latitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_den_col(), dtype=np.float64)
    )
    res_py_den_row = derivative_polynomial_latitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_den_row(), dtype=np.float64)
    )
    res_py_num_col = derivative_polynomial_latitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_num_col(), dtype=np.float64)
    )
    res_py_num_row = derivative_polynomial_latitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_num_row(), dtype=np.float64)
    )

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
def test_derivative_polynomial_longitude(geom_path):
    """
    test derivative_polynomial_longitude function
    """

    file_path = os.path.join(data_path(), geom_path)
    rpc_cpp = bindings_cpp_constructor(file_path)

    # arbitrary values (extract from a rpc.py test)
    xnorm = -0.95821893  # lon_norm
    ynorm = 0.97766941  # lat_norm
    znorm = -5.29411765  # alt_norm

    res_c_den_col = bindings_cpp.derivative_polynomial_longitude(xnorm, ynorm, znorm, rpc_cpp.get_den_col())
    res_c_den_row = bindings_cpp.derivative_polynomial_longitude(xnorm, ynorm, znorm, rpc_cpp.get_den_row())
    res_c_num_col = bindings_cpp.derivative_polynomial_longitude(xnorm, ynorm, znorm, rpc_cpp.get_num_col())
    res_c_num_row = bindings_cpp.derivative_polynomial_longitude(xnorm, ynorm, znorm, rpc_cpp.get_num_row())

    # rpc.py derivative_polynomial_longitude
    res_py_den_col = derivative_polynomial_longitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_den_col(), dtype=np.float64)
    )
    res_py_den_row = derivative_polynomial_longitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_den_row(), dtype=np.float64)
    )
    res_py_num_col = derivative_polynomial_longitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_num_col(), dtype=np.float64)
    )
    res_py_num_row = derivative_polynomial_longitude(
        xnorm, ynorm, znorm, np.array(rpc_cpp.get_num_row(), dtype=np.float64)
    )

    assert res_c_den_col == pytest.approx(res_py_den_col, abs=1e-15)
    assert res_c_den_row == pytest.approx(res_py_den_row, abs=1e-15)
    assert res_c_num_col == pytest.approx(res_py_num_col, abs=1e-15)
    assert res_c_num_row == pytest.approx(res_py_num_row, abs=1e-15)


def test_compute_loc_inverse_derivates():
    """
    test on the compute_loc_inverse_derivates methode
    """

    rpc_path = os.path.join(data_path(), "rpc/PHRDIMAP_P1BP--2018122638935449CP.XML")

    rpc_optim = GeoModel(rpc_path, "RPCoptim")
    rpc_py = GeoModel(rpc_path, "RPC")

    # INPUTS
    nrb_point = 1e3
    first_lon = 7.0477886581984
    first_lat = 43.62208491280199
    last_lon = 7.308411551163017
    last_lat = 43.73298365695963
    first_alt = -50
    last_alt = 1000

    lon_vect = np.linspace(first_lon, last_lon, int(nrb_point ** (1 / 3)) + 1)
    lat_vect = np.linspace(first_lat, last_lat, int(nrb_point ** (1 / 3)) + 1)
    alt_vect = np.linspace(first_alt, last_alt, int(nrb_point ** (1 / 3)) + 1)

    lon_vect, lat_vect, alt_vect = np.meshgrid(lon_vect, lat_vect, alt_vect)

    lon_vect = np.ndarray.flatten(lon_vect)
    lat_vect = np.ndarray.flatten(lat_vect)
    alt_vect = np.ndarray.flatten(alt_vect)

    res_py = rpc_py.compute_loc_inverse_derivates(lon_vect, lat_vect, alt_vect)
    res_cpp = np.empty((len(lon_vect), 4))
    for i, lon_vect_i in enumerate(lon_vect):
        res_cpp[i, :] = rpc_optim.compute_loc_inverse_derivates(lon_vect_i, lat_vect[i], alt_vect[i])

    np.testing.assert_allclose(res_cpp[:, 0], res_py[0], 0, 6e-10)
    np.testing.assert_allclose(res_cpp[:, 1], res_py[1], 0, 6e-10)
    np.testing.assert_allclose(res_cpp[:, 2], res_py[2], 0, 2e-10)
    np.testing.assert_allclose(res_cpp[:, 3], res_py[3], 0, 3e-10)


@pytest.mark.parametrize(
    "col,row,alt",
    [
        (600.0, 200.0, 125.0),
        (0.0, 0.0, 0.0),
        (150.0, 300.0, 50.0),
        (354.0, 124.0, 101.0),
        (537.3, 178.7, 300.8),
    ],
)
def test_rpc_direct_inverse_iterative_unitary_loc(col, row, alt):
    """
    test direct_loc_inverse_iterative methode on single localisation
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    rpc_cpp = GeoModel(file_dimap, "RPCoptim")
    rpc_py = GeoModel(file_dimap, "RPC")

    nb_iter_max = 10
    fill_nan = False

    (lon_py, lat_py, alt_py) = rpc_py.direct_loc_inverse_iterative(row, col, alt, nb_iter_max, fill_nan)
    (lon_cpp, lat_cpp, alt_cpp) = rpc_cpp.direct_loc_inverse_iterative([row], [col], [alt], nb_iter_max, fill_nan)

    assert lon_cpp == lon_py[0]
    assert lat_cpp == lat_py[0]
    assert alt_cpp == alt_py[0]

    # Check identity

    (row_cpp, col_cpp, alt_cpp) = rpc_cpp.inverse_loc(lon_cpp, lat_cpp, alt_cpp)

    # small error from numba
    assert row_cpp[0] == pytest.approx(row, abs=1e-9)
    assert col_cpp[0] == pytest.approx(col, abs=1e-9)
    assert alt_cpp[0] == pytest.approx(alt, abs=0)


def test_rpc_direct_inverse_iterative_multi_loc():
    """
    test direct_loc_inverse_iterative methode on multiple localisations
    """
    rpc_path = os.path.join(data_path(), "rpc/PHRDIMAP_P1BP--2018122638935449CP.XML")

    rpc_cpp = GeoModel(rpc_path, "RPCoptim")
    rpc_py = GeoModel(rpc_path, "RPC")

    # INPUTS
    nrb_point = 1e4
    first_row = 1
    first_col = 1
    last_row = 22940
    last_col = 40000

    row_vect = np.linspace(first_row, last_row, int(nrb_point ** (1 / 2)))
    col_vect = np.linspace(first_col, last_col, int(nrb_point ** (1 / 2)))

    row_vect, col_vect = np.meshgrid(row_vect, col_vect)

    row_vect = np.ndarray.flatten(row_vect)
    col_vect = np.ndarray.flatten(col_vect)
    alt_vect = np.linspace(-100, 10000, len(col_vect))

    nb_iter_max = 10
    fill_nan = False

    (lon_py_vect, lat_py_vect, alt_py_vect) = rpc_py.direct_loc_inverse_iterative(
        row_vect, col_vect, alt_vect, nb_iter_max, fill_nan
    )
    (lon_cpp_vect, lat_cpp_vect, alt_cpp_vect) = rpc_cpp.direct_loc_inverse_iterative(
        row_vect, col_vect, alt_vect, nb_iter_max, fill_nan
    )

    np.testing.assert_allclose(np.array(lon_cpp_vect), lon_py_vect, 0, 1e-11)
    np.testing.assert_allclose(np.array(lat_cpp_vect), lat_py_vect, 0, 1e-11)
    np.testing.assert_allclose(np.array(alt_cpp_vect), alt_py_vect, 0, 1e-11)

    # --- Differents lenght ---#
    row_vect_min = row_vect[int(len(row_vect) / 2) :]
    (lon_py, lat_py, alt_py) = rpc_py.direct_loc_inverse_iterative(
        row_vect_min, col_vect[: int(len(col_vect) / 2)], alt_vect, nb_iter_max, fill_nan
    )
    (lon_cpp, lat_cpp, alt_cpp) = rpc_cpp.direct_loc_inverse_iterative(
        row_vect_min, col_vect, alt_vect, nb_iter_max, fill_nan
    )

    np.testing.assert_allclose(np.array(lon_cpp), lon_py, 0, 1e-11)
    np.testing.assert_allclose(np.array(lat_cpp), lat_py, 0, 1e-11)
    np.testing.assert_allclose(np.array(alt_cpp), alt_py, 0, 1e-11)

    # --- Nan ---#
    row_vect_nan = row_vect
    row_vect_nan[1::2] = np.nan

    (lon_py, lat_py, alt_py) = rpc_py.direct_loc_inverse_iterative(
        row_vect_nan, col_vect, alt_vect, nb_iter_max, fill_nan
    )
    (lon_cpp, lat_cpp, alt_cpp) = rpc_cpp.direct_loc_inverse_iterative(
        row_vect_nan, col_vect, alt_vect, nb_iter_max, fill_nan
    )

    np.testing.assert_allclose(np.array(lon_cpp), lon_py, 0, 1e-11)
    np.testing.assert_allclose(np.array(lat_cpp), lat_py, 0, 1e-11)
    np.testing.assert_allclose(np.array(alt_cpp), alt_py, 0, 1e-11)

    # --- Full Nan ---#
    nans = np.full((10), np.nan, dtype=float)
    res_py = rpc_py.direct_loc_inverse_iterative(nans, nans, nans, nb_iter_max, fill_nan)
    res_cpp = rpc_cpp.direct_loc_inverse_iterative(nans, nans, nans, nb_iter_max, fill_nan)

    res_py = np.vstack([res_py[0], res_py[1]])
    res_cpp = np.array(res_cpp[0:2])

    if np.isnan(res_py).any() and np.isnan(res_cpp).any():
        assert True
    else:
        assert AssertionError()


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
def test_direct_loc_h(geom_path):
    """
    test test_direct_loc_h
    """

    rpc_path = os.path.join(data_path(), geom_path)

    rpc_optim = GeoModel(rpc_path, "RPCoptim")
    rpc_py = GeoModel(rpc_path, "RPC")

    # INPUTS
    nrb_point = 1e6
    first_row = 1
    first_col = 1
    last_row = 22940
    last_col = 40000
    first_alt = -50
    last_alt = 1000

    col_vect = np.linspace(first_col, last_col, int(nrb_point ** (1 / 3)) + 1)
    row_vect = np.linspace(first_row, last_row, int(nrb_point ** (1 / 3)) + 1)
    alt_vect = np.linspace(first_alt, last_alt, int(nrb_point ** (1 / 3)) + 1)

    col_vect, row_vect, alt_vect = np.meshgrid(col_vect, row_vect, alt_vect)

    col_vect = np.ndarray.flatten(col_vect)
    row_vect = np.ndarray.flatten(row_vect)
    alt_vect = np.ndarray.flatten(alt_vect)

    col_vect[0] = np.nan
    res_optim = rpc_optim.direct_loc_h(row_vect[:5], col_vect[:5], alt_vect[:5])
    res_py = rpc_py.direct_loc_h(row_vect[:5], col_vect[:5], alt_vect[:5])

    np.testing.assert_allclose(res_optim[:, 0], res_py[:, 0], 0, 9e-16)
    np.testing.assert_allclose(res_optim[:, 1], res_py[:, 1], 0, 0)
    np.testing.assert_allclose(res_optim[:, 2], res_py[:, 2], 0, 0)


def test_get_alt_min_max():
    """ "
    test get_alt_min_max method
    """

    rpc_path = os.path.join(data_path(), "rpc/PHRDIMAP_P1BP--2018122638935449CP.XML")

    rpc_cpp = bindings_cpp_constructor(rpc_path)
    rpc_py = GeoModel(rpc_path, "RPC")

    assert rpc_py.get_alt_min_max() == rpc_cpp.get_alt_min_max()


def test_los_extrema():
    """
    test los_extrema method
    """

    rpc_path = os.path.join(data_path(), "rpc/PHRDIMAP_P1BP--2018122638935449CP.XML")

    rpc_cpp = bindings_cpp_constructor(rpc_path)
    rpc_py = GeoModel(rpc_path, "RPC")
    rpc_optim = GeoModel(rpc_path, "RPCoptim")

    row = 100.5
    col = 200.5
    alt_min = -10
    alt_max = 2000
    fill_nan = False
    epsg = 4326

    # --Normal case
    res_py = rpc_py.los_extrema(row, col, alt_min, alt_max, fill_nan, epsg)
    res_cpp = rpc_cpp.los_extrema(row, col, alt_min, alt_max, fill_nan)
    res_optim = rpc_optim.los_extrema(row, col, alt_min, alt_max, fill_nan)

    res_cpp = np.array(res_cpp).T

    np.testing.assert_allclose(res_py, res_cpp, 0, 0)
    np.testing.assert_allclose(res_py, res_optim, 0, 0)

    # --Case alt_max=None/nan for c++
    res_py = rpc_py.los_extrema(row, col, alt_min, None, fill_nan, epsg)
    res_cpp = rpc_cpp.los_extrema(row, col, alt_min, np.nan, fill_nan)
    res_optim = rpc_optim.los_extrema(row, col, alt_min, np.nan, fill_nan)

    res_cpp = np.array(res_cpp).T

    np.testing.assert_allclose(res_py, res_cpp, 0, 0)
    np.testing.assert_allclose(res_py, res_optim, 0, 0)

    # --Case alt_min >= alt_minmax[0] && alt_max <= alt_minmax[1]
    alt_min = 15000
    alt_max = -15000
    res_py = rpc_py.los_extrema(row, col, alt_min, alt_max, fill_nan, epsg)
    res_cpp = rpc_cpp.los_extrema(row, col, alt_min, alt_max, fill_nan)
    res_optim = rpc_optim.los_extrema(row, col, alt_min, alt_max, fill_nan)

    res_cpp = np.array(res_cpp).T

    np.testing.assert_allclose(res_py, res_cpp, 0, 0)
    np.testing.assert_allclose(res_py, res_optim, 0, 0)

    # -- case epsg!=4326
    epsg = 32632
    res_py = rpc_py.los_extrema(row, col, alt_min, alt_max, fill_nan, epsg)
    res_optim = rpc_optim.los_extrema(row, col, alt_min, alt_max, fill_nan, epsg)

    np.testing.assert_allclose(res_py, res_optim, 0, 0)


def test_direct_loc_dtm():
    """
    test direct_loc_dtm method
    """

    rpc_file = os.path.join(data_path(), "rpc/phr_ventoux/RPC_PHR1B_P_201308051042194_SEN_690908101-001.XML")
    mnt = os.path.join(data_path(), "dtm/srtm_ventoux/srtm90_non_void_filled/N44E005.hgt")

    rpc_py = GeoModel(rpc_file, "RPC")
    rpc_optim = GeoModel(rpc_file, "RPCoptim")

    dtm_image = dtm_reader(mnt)
    dtm_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    dtm_cpp = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    nrb_point = 1e4
    first_row = 1
    first_col = 1
    last_row = 38608
    last_col = 36416

    row_vect = np.linspace(first_row, last_row, int(nrb_point ** (1 / 2)) + 1)
    col_vect = np.linspace(first_col, last_col, int(nrb_point ** (1 / 2)) + 1)

    row_vect, col_vect = np.meshgrid(row_vect, col_vect)

    row_vect = np.ndarray.flatten(row_vect)
    col_vect = np.ndarray.flatten(col_vect)

    col_vect[0] = np.nan

    res_py = rpc_py.direct_loc_dtm(row_vect, col_vect, dtm_py)
    res_optim = rpc_optim.direct_loc_dtm(row_vect, col_vect, dtm_cpp)

    np.testing.assert_allclose(res_optim[:, 0], res_py[:, 0], 0, 9e-12)
    np.testing.assert_allclose(res_optim[:, 1], res_py[:, 1], 0, 6e-12)
    np.testing.assert_allclose(res_optim[:, 2], res_py[:, 2], 0, 5e-7)


def test_direct_loc_dtm_utm():
    """
    Test direct loc dtm
    """

    rpc_file = os.path.join(data_path(), "rpc/phr_ventoux/RPC_PHR1B_P_201308051042194_SEN_690908101-001.XML")
    mnt = os.path.join(data_path(), "dtm/srtm_ventoux/srtm90_resampled_UTM31/N44E005_UTM.tif")

    rpc_py = GeoModel(rpc_file, "RPC")
    rpc_optim = GeoModel(rpc_file, "RPCoptim")

    dtm_image = dtm_reader(mnt)
    dtm_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    dtm_cpp = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    nrb_point = 1e4
    first_row = 1
    first_col = 1
    last_row = 38608
    last_col = 36416

    row_vect = np.linspace(first_row, last_row, int(nrb_point ** (1 / 2)) + 1)
    col_vect = np.linspace(first_col, last_col, int(nrb_point ** (1 / 2)) + 1)

    row_vect, col_vect = np.meshgrid(row_vect, col_vect)

    row_vect = np.ndarray.flatten(row_vect)
    col_vect = np.ndarray.flatten(col_vect)

    res_py = rpc_py.direct_loc_dtm(row_vect, col_vect, dtm_py)
    res_optim = rpc_optim.direct_loc_dtm(row_vect, col_vect, dtm_cpp)
    res_py_dtm_optim = rpc_py.direct_loc_dtm(row_vect, col_vect, dtm_cpp)

    np.testing.assert_array_equal(res_py, res_py_dtm_optim)
    np.testing.assert_allclose(res_optim[:, 0], res_py[:, 0], 0, 7e-7)
    np.testing.assert_allclose(res_optim[:, 1], res_py[:, 1], 0, 7e-7)
    np.testing.assert_allclose(res_optim[:, 2], res_py[:, 2], 0, 4e-7)
