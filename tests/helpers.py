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
Helpers shared testing generic module:
contains global shared generic functions for tests/*.py
"""
# pylint: disable=R0801, C0103

# Standard imports
import os

# Shareloc bindings
import bindings_cpp
from shareloc.dtm_reader import dtm_reader
from shareloc.geofunctions.dtm_intersection import DTMIntersection

# Shareloc imports
from shareloc.geomodels.rpc_readers import rpc_reader


def data_path(alti="", scene_id=""):
    """
    return the data path, when used without any argument data_path() returns data directory
    :param alti: first sub dir corresponding to datum ("ellipsoide" or "geoid")
    :type alti: str
    :param scene_id: second sub dir corresponding to the scene id
    :type scene_id: str
    :return: data path.
    :rtype: str
    """
    data_root_folder = os.path.join(os.path.dirname(__file__), "data")
    sub_folder = os.path.join(alti, scene_id)
    return os.path.join(data_root_folder, sub_folder)


def bindings_cpp_constructor(file_path: str):
    """
    Create rpc c++ object from rpc file path

    :param file_path: path to rpc file
    :type file_path: str
    :return: RPC c++ object
    :retype: bindings_cpp.RPC
    """

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

    empty = [0 for i in range(20)]

    if rpc_params["num_col"] and rpc_params["num_x"]:  # direct and inverse coef
        rpc_cpp = bindings_cpp.RPC(
            True,
            True,
            rpc_params["num_col"],
            rpc_params["den_col"],
            rpc_params["num_row"],
            rpc_params["den_row"],
            rpc_params["num_x"],
            rpc_params["den_x"],
            rpc_params["num_y"],
            rpc_params["den_y"],
            norm_coeffs,
        )
    elif rpc_params["num_col"]:  # only inverse coef
        rpc_cpp = bindings_cpp.RPC(
            True,
            False,
            rpc_params["num_col"],
            rpc_params["den_col"],
            rpc_params["num_row"],
            rpc_params["den_row"],
            empty,
            empty,
            empty,
            empty,
            norm_coeffs,
        )
    elif rpc_params["num_x"]:  # only direct coef
        rpc_cpp = bindings_cpp.RPC(
            False,
            True,
            empty,
            empty,
            empty,
            empty,
            rpc_params["num_x"],
            rpc_params["den_x"],
            rpc_params["num_y"],
            rpc_params["den_y"],
            norm_coeffs,
        )
    else:
        raise ValueError("RPCoptim : No RPC coefficients readable")

    return rpc_cpp


def DTMIntersection_constructor(file_path: str):
    """
    Create DTMIntersection python and c++ object from mnt file path

    :param file_path: path to rpc file
    :type file_path: str
    :return: python and c++ DTMIntersection object
    :retype: tuple(DTMIntersection python, DTMIntersection c++)
    """
    dtm_image = dtm_reader(file_path)
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
    return dtm_py, dtm_cpp
