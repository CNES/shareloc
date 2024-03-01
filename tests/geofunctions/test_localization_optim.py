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
Test module for localisation class shareloc/geofunctions/localisation_optim.py
"""
# Standard imports
import os

# Third party imports
import numpy as np

# Shareloc bindings
import bindings_cpp
from shareloc.geofunctions.localization import coloc
from shareloc.geomodels import GeoModel

# Shareloc test imports
from ..helpers import DTMIntersection_constructor, bindings_cpp_constructor, data_path


def test_coloc():
    """
    Test colocation function
    """
    mnt = os.path.join(data_path(), "dtm/srtm_ventoux/srtm90_non_void_filled/N44E005.hgt")

    data_left = os.path.join(data_path(), "rectification", "left_image")
    rpc_py_left = GeoModel(data_left + ".geom", "RPC")
    rpc_cpp_left = bindings_cpp_constructor(data_left + ".geom")

    data_right = os.path.join(data_path(), "rectification", "right_image")
    rpc_py_right = GeoModel(data_right + ".geom", "RPC")
    rpc_cpp_right = bindings_cpp_constructor(data_right + ".geom")

    dtm_py, dtm_cpp = DTMIntersection_constructor(mnt)

    row_vect = np.array([100, 200])
    col_vect = np.array([100, 200])
    alt_vect = np.array([1275, 1475])

    # Identity

    res_id_1 = bindings_cpp.coloc(rpc_cpp_left, rpc_cpp_left, row_vect, col_vect, dtm_cpp)
    res_id_2 = bindings_cpp.coloc(rpc_cpp_left, rpc_cpp_left, row_vect, col_vect, alt_vect)

    np.testing.assert_allclose(row_vect, res_id_1[0], 0, 3e-2)
    np.testing.assert_allclose(col_vect, res_id_1[1], 0, 8e-3)
    np.testing.assert_allclose(row_vect, res_id_2[0], 0, 6e-10)
    np.testing.assert_allclose(col_vect, res_id_2[1], 0, 5e-11)

    # compare to python

    res_cpp_1 = bindings_cpp.coloc(rpc_cpp_right, rpc_cpp_left, row_vect, col_vect, dtm_cpp)
    res_cpp_2 = bindings_cpp.coloc(rpc_cpp_right, rpc_cpp_left, row_vect, col_vect, alt_vect)

    res_py_1 = coloc(rpc_py_right, rpc_py_left, row_vect, col_vect, dtm_py, using_geotransform=False)
    res_py_2 = coloc(rpc_py_right, rpc_py_left, row_vect, col_vect, alt_vect, using_geotransform=False)

    np.testing.assert_allclose(res_cpp_1[0], res_py_1[0], 0, 2e-11)
    np.testing.assert_allclose(res_cpp_1[1], res_py_1[1], 0, 4e-12)
    np.testing.assert_allclose(res_cpp_1[2], res_py_1[2], 0, 0)

    np.testing.assert_allclose(res_cpp_2[0], res_py_2[0], 0, 2e-11)
    np.testing.assert_allclose(res_cpp_2[1], res_py_2[1], 0, 2e-11)
    np.testing.assert_allclose(res_cpp_2[2], res_py_2[2], 0, 0)
