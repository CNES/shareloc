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
Module to test geomodel_template class loading capabilities
"""
# TODO: refactor to disable no-member in RPC class
# pylint: disable=no-member


# Standard imports
import os

# Third party imports
import pytest
import rasterio
import json

# Shareloc imports
from shareloc.geomodels.geomodel_template import GeoModelTemplate, identify_dimap, identify_ossim_kwl

# Shareloc test imports
from ..helpers import data_path


def test_rpc_drivers():
    """
    test GeoModelTemplate driver identification
    """
    data_folder = data_path()

    # Test DIMAP GeoModelTemplate
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")
    fctrat_dimap = GeoModelTemplate(file_dimap)
    assert fctrat_dimap.rpc_params["driver_type"] == "dimap_v1.4"

    # Test OSSIM KWL GEOM GeoModelTemplate
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    fctrat_geom = GeoModelTemplate(file_geom)
    assert fctrat_geom.rpc_params["driver_type"] == "ossim_kwl"

    # Test DIMAPv2 GeoModelTemplate
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_dimap_v2 = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    fctrat_dimap_v2 = GeoModelTemplate(file_dimap_v2)
    assert fctrat_dimap_v2.rpc_params["driver_type"] == "dimap_v2.15"

    # Test fake GeoModelTemplate
    fake_rpc = os.path.join(data_folder, "rpc/fake_rpc.txt")
    try:
        GeoModelTemplate(fake_rpc)  # Raise ValueError->True
        raise AssertionError()  # Assert false if no exception raised
    except ValueError:
        assert True


def test_identify_ossim_kwl():
    """
    test ossim file identification
    """
    data_folder = data_path()
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    ossim_model = identify_ossim_kwl(file_geom)
    assert ossim_model == "ossimPleiadesModel"


def test_identify_dimap():
    """
    test dimap file identification
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")
    dimap_version = identify_dimap(file_dimap)
    assert dimap_version == "1.4"



@pytest.mark.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)
@pytest.mark.parametrize(
    "prod, can_read",
    [
        ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif", True),
        ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001_nogeo.tif", False),
    ],
)
def test_GeoModelTemplate_from_geotiff_without_rpc(prod, can_read):
    """
    test geotiff file without rpc
    """
    data_folder = data_path()
    rpc_file = os.path.join(data_folder, "rpc", prod)
    try:
        GeoModelTemplate(rpc_file)
        assert can_read
    except ValueError:
        assert not can_read


def test_load_rpc_params():
    """
    test loading of rpc_params dict
    """
    geom = GeoModelTemplate("tests/data/rpc/PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.geom").rpc_params
    tif = GeoModelTemplate("tests/data/rpc/PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif").rpc_params
    PHRDIMAP = GeoModelTemplate("tests/data/rpc/PHRDIMAP_P1BP--2017030824934340CP.XML").rpc_params
    RPC_P1BP = GeoModelTemplate("tests/data/rpc/RPC_P1BP--2017092838284574CP.XML").rpc_params
    RPC_PHR1B = GeoModelTemplate("tests/data/rpc/RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML").rpc_params


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