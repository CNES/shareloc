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
# pylint: disable=no-member, invalid-name


# Third party imports
import json

from shareloc.geomodels import GeoModel
from shareloc.geomodels.rpc_optim import RpcOptim  # noqa : F401 # pylint: disable=unused-import


def test_load_rpc_params():
    """
    test loading of rpc_params dict
    """
    geom = GeoModel("tests/data/rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom", "RpcOptim").__dict__
    tif = GeoModel("tests/data/rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif", "RpcOptim").__dict__
    PHRDIMAP = GeoModel("tests/data/rpc/PHRDIMAP_P1BP--2017030824934340CP.XML", "RpcOptim").__dict__
    RPC_P1BP = GeoModel("tests/data/rpc/RPC_P1BP--2017092838284574CP.XML", "RpcOptim").__dict__
    RPC_PHR1B = GeoModel("tests/data/rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML", "RpcOptim").__dict__

    del geom["type"]
    del tif["type"]
    del PHRDIMAP["type"]
    del RPC_P1BP["type"]
    del RPC_PHR1B["type"]

    with open("tests/data/geomodel_template/geom.json", "r", encoding="utf-8") as f:
        geom_ref = json.load(f)
    with open("tests/data/geomodel_template/tif.json", "r", encoding="utf-8") as f:
        tif_ref = json.load(f)
    with open("tests/data/geomodel_template/PHRDIMAP.json", "r", encoding="utf-8") as f:
        PHRDIMAP_ref = json.load(f)
    with open("tests/data/geomodel_template/RPC_P1BP.json", "r", encoding="utf-8") as f:
        RPC_P1BP_ref = json.load(f)
    with open("tests/data/geomodel_template/RPC_PHR1B.json", "r", encoding="utf-8") as f:
        RPC_PHR1B_ref = json.load(f)

    assert geom == geom_ref
    assert tif == tif_ref
    assert PHRDIMAP == PHRDIMAP_ref
    assert RPC_P1BP == RPC_P1BP_ref
    assert RPC_PHR1B == RPC_PHR1B_ref


def test_function_rpc_cpp():
    """
    Test call to methode parent class rpc in cpp
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
