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
Module to test RPC_optim class
"""
# pylint: disable=no-member


# Third party imports
import json

# Shareloc imports
from shareloc.geomodels.rpc_optim import RPC_optim




def test_load_rpc_params():
    """
    test loading of rpc_params dict
    """
    geom = RPC_optim("tests/data/rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom").rpc_params
    tif = RPC_optim("tests/data/rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif").rpc_params
    PHRDIMAP = RPC_optim("tests/data/rpc/PHRDIMAP_P1BP--2017030824934340CP.XML").rpc_params
    RPC_P1BP = RPC_optim("tests/data/rpc/RPC_P1BP--2017092838284574CP.XML").rpc_params
    RPC_PHR1B = RPC_optim("tests/data/rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML").rpc_params


    a = open("tests/data/geomodel_template/geom.json")
    geom_ref = json.load(a)
    b = open("tests/data/geomodel_template/tif.json")
    tif_ref = json.load(b)
    c = open("tests/data/geomodel_template/PHRDIMAP.json")
    PHRDIMAP_ref = json.load(c)
    d = open("tests/data/geomodel_template/RPC_P1BP.json")
    RPC_P1BP_ref = json.load(d)
    e = open("tests/data/geomodel_template/RPC_PHR1B.json")
    RPC_PHR1B_ref = json.load(e)

    assert geom == geom_ref
    assert tif == tif_ref
    assert PHRDIMAP == PHRDIMAP_ref
    assert RPC_P1BP == RPC_P1BP_ref
    assert RPC_PHR1B == RPC_PHR1B_ref

    a.close()
    b.close()
    c.close()
    d.close()
    e.close()

def test_function_rpc_cpp():
    """
    Test call to methode parent class rpc in cpp
    """

    rpc = RPC_optim("tests/data/rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom")

    vector_double = [1.,1.,1.]
    double = 1.
    int = 1
    string = "peuimporte"

    try:
        rpc.direct_loc_h(vector_double,vector_double, double, False)
        rpc.direct_loc_grid_h(int, int, int, int, int, int, double )
        rpc.direct_loc_dtm(double, double, string)
        rpc.inverse_loc(vector_double, vector_double, double)
        rpc.filter_coordinates(vector_double, vector_double, False, string )
        rpc.compute_loc_inverse_derivates(vector_double, vector_double, vector_double)
        rpc.direct_loc_inverse_iterative(vector_double, vector_double, double, int, False)
        rpc.get_alt_min_max()
        rpc.los_extrema(double, double, double , double, False, int)

    except:
        assert False

    assert True