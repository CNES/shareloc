#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
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
Module to test functions that use rpc
"""


import os
import pytest
from utils import test_path
import numpy as np

from shareloc.rpc.rpc import RPC, identify_dimap, identify_ossim_kwl, identify_euclidium_rpc


def test_rpc_drivers():
    """
    test rpc driver identification
    """
    data_folder = test_path()

    file_direct_euclide = os.path.join(data_folder, "rpc/rpc_dir.euc")
    file_inverse_euclide = os.path.join(data_folder, "rpc/rpc_inv.euc")

    fctrat_eucl = RPC.from_any(file_inverse_euclide, file_direct_euclide)

    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat_dimap = RPC.from_any(file_dimap)

    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, "rpc/{}.geom".format(id_scene))
    fctrat_geom = RPC.from_any(file_geom)

    file_dimap = os.path.join(data_folder, "rpc/RPC_{}.XML".format(id_scene))
    print(file_dimap)
    fctrat_dimap_v2 = RPC.from_any(file_dimap, topleftconvention=True)

    fake_rpc = os.path.join(data_folder, "rpc/fake_rpc.txt")
    fctrat_fake = RPC.from_any(fake_rpc, topleftconvention=True)
    assert fctrat_eucl.driver_type == "euclidium"
    assert fctrat_dimap.driver_type == "dimap_v1.4"
    assert fctrat_geom.driver_type == "ossim_kwl"
    assert fctrat_dimap_v2.driver_type == "dimap_v2.15"
    assert fctrat_fake is None


def test_identify_ossim_kwl():
    """
    test ossim file identification
    """
    data_folder = test_path()
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, "rpc/{}.geom".format(id_scene))
    ossim_model = identify_ossim_kwl(file_geom)
    assert ossim_model == "ossimPleiadesModel"


def test_identify_euclidium_rpc():
    """
    test euclidium file identification
    """
    data_folder = test_path()
    file_direct_euclide = os.path.join(data_folder, "rpc/rpc_dir.euc")
    file_inverse_euclide = os.path.join(data_folder, "rpc/rpc_inv.euc")

    is_eucl_inverse = identify_euclidium_rpc(file_inverse_euclide)
    is_eucl_direct = identify_euclidium_rpc(file_direct_euclide)

    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, "rpc/{}.geom".format(id_scene))
    is_eucl_bad = identify_euclidium_rpc(file_geom)
    assert is_eucl_inverse
    assert is_eucl_direct
    assert not is_eucl_bad


def test_identify_dimap():
    """
    test dimap file identification
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))
    dimap_version = identify_dimap(file_dimap)
    assert dimap_version == "1.4"


@pytest.mark.parametrize("id_scene,lon,lat,alt, col_vt,row_vt", [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001",
                                                                  7.048662660737769592, 43.72774839443545858, 0.0,
                                                                  100.0, 200.0),
                                                                 ("PHR1B_P_201709281038045_SEN_PRG_FC_178608-001",
                                                                  7.11526088296757386331137240632,
                                                                  43.684281179313565246502548689, 850.0,  10121.0657,
                                                                  10235.9310)])
def test_rpc_ossim_kwl(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization from ossim file
    """
    data_folder = test_path()
    file_geom = os.path.join(data_folder, "rpc/{}.geom".format(id_scene))
    fctrat_geom = RPC.from_any(file_geom, topleftconvention=True)
    print("{} {} {}".format(lon, lat, alt))
    (row, col, __) = fctrat_geom.inverse_loc(lon, lat, alt)
    print("col {} row {}".format(col, row))
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)


@pytest.mark.parametrize("id_scene,lon,lat,alt, col_vt,row_vt",
                         [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001_inv.txt", 7.048662660737769592,
                           43.72774839443545858, 0.0, 100.0, 200.0),
                          ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.geom", 7.048662660737769592,
                           43.72774839443545858, 0.0, 100.0, 200.0),
                          ("RPC_PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.XML", 7.048662660737769592,
                           43.72774839443545858, 0.0, 100.0, 200.0),
                          ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif", 7.048662660737769592,
                           43.72774839443545858, 0.0, 100.0, 200.0), ])
def test_rpc_from_any(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization from any file
    """
    data_folder = test_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)
    fctrat = RPC.from_any(rpc_file, topleftconvention=True)
    print("{} {} {}".format(lon, lat, alt))
    (row, col, __) = fctrat.inverse_loc(lon, lat, alt)
    print("col {} row {}".format(col, row))
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)


@pytest.mark.parametrize("id_scene,lon,lat,alt", [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001_inv.txt",
                                                   7.048662660737769592, 43.72774839443545858, 0.0)])
def test_rpc_euclidium_direct_iterative(id_scene, lon, lat, alt):
    """
    test the sequence of a inverse localization followed by a direct localization using euclidium file
    """
    data_folder = test_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)
    fctrat = RPC.from_any(rpc_file, topleftconvention=True)
    print("{} {} {}".format(lon, lat, alt))
    (row, col, __) = fctrat.inverse_loc(lon, lat, alt)
    print("col {} row {}".format(col, row))
    (lon2, lat2, __) = fctrat.direct_loc_h(row, col, alt)
    assert lon == pytest.approx(lon2, abs=1e-2)
    assert lat == pytest.approx(lat2, abs=1e-2)


@pytest.mark.parametrize("id_scene,lon,lat,alt, col_vt,row_vt",
                         [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001_dir.txt", 7.048662660737769592,
                           43.72774839443545858, 0.0, 100.0, 200.0)])
def test_rpc_euclidium_direct(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test direct and inverse localization using fake rpc file
    """
    data_folder = test_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)
    fctrat = RPC.from_any(rpc_file, topleftconvention=True)
    (row, col, __) = fctrat.inverse_loc(lon, lat, alt)
    assert row is None
    assert col is None
    (lon2, lat2, __) = fctrat.direct_loc_inverse_iterative(row_vt, col_vt, alt)
    assert lat2 is None
    assert lon2 is None


@pytest.mark.parametrize("prod, can_read", [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif", True),
                                            ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001_nogeo.tif", False)])
def test_rpc_from_geotiff_without_rpc(prod, can_read):
    """
    test file without rpc
    """
    data_folder = test_path()
    rpc_file = os.path.join(data_folder, "rpc",  prod)
    try:
        RPC.from_geotiff(rpc_file, topleftconvention=True)
        assert can_read
    except ValueError:
        assert not can_read


@pytest.mark.parametrize("id_scene,lon,lat,alt, col_vt,row_vt", [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001",
                                                                  7.048662660737769592, 43.72774839443545858, 0.0,
                                                                  100.0, 200.0)])
def test_rpc_dimap_v2(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization using ossim and dimap files
    """
    data_folder = test_path()
    file_dimap = os.path.join(data_folder, "rpc/RPC_{}.XML".format(id_scene))
    print(file_dimap)
    fctrat_dimap = RPC.from_dimap(file_dimap, topleftconvention=True)
    print("{} {} {}".format(lon, lat, alt))
    (row, col, __) = fctrat_dimap.inverse_loc(lon, lat, alt)
    print("col {} row {} dimap".format(col, row))
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)

    file_geom = os.path.join(data_folder, "rpc/{}.geom".format(id_scene))
    fctrat_geom = RPC.from_any(file_geom, topleftconvention=True)
    (row, col, __) = fctrat_geom.inverse_loc(lon, lat, alt)
    print("col {} row {} geom".format(col, row))
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_phrdimap(col, row, alt):
    """
    test the sequence of a inverse localization followed by a direct localization using dimap file
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_dimap(file_dimap)

    (lon, lat, alt) = fctrat.direct_loc_h(row, col, alt)
    print (lon, lat)

    (row_ar, col_ar, __) = fctrat.inverse_loc(lon, lat, alt)
    print(col_ar, row_ar)
    assert col_ar == pytest.approx(col, abs=2e-2)
    assert row_ar == pytest.approx(row, abs=2e-2)


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_direct_inverse_iterative_vs_direct(col, row, alt):
    """
    compare direct localization and iterative direct localization
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    #(col,lig,alt)=(100,1000,400)
    (x_direct, y_direct, __) = fctrat.direct_loc_h(row, col, alt)
    (x_inv, y_inv, __) = fctrat.direct_loc_inverse_iterative(row, col, alt)
    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    assert x_direct == pytest.approx(x_inv, abs= 10.0/111111000)
    assert y_direct == pytest.approx(y_inv, abs= 10.0/111111000)


def test_rpc_direct_inverse_iterative_vs_direct_multiple_points():
    """
    compare direct localization and iterative direct localization with multiple points
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (col, row, alt) = (np.array([600, 610]), np.array([200, 210]), np.array([125]))
    p_direct = fctrat.direct_loc_h(row, col, alt)
    (rowinv, colinv, __) = fctrat.inverse_loc(p_direct[:, 0], p_direct[:, 1], alt)
    p_direct_iterative = fctrat.direct_loc_inverse_iterative(rowinv, colinv, alt)

    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    #Point error col = 600, row = 200
    assert p_direct[0, 0] == p_direct_iterative[0][0]
    assert p_direct[0, 1] == p_direct_iterative[1][0]

    # Point error col = 601, row = 201
    assert p_direct[1, 0] == p_direct_iterative[0][1]
    assert p_direct[1, 1] == p_direct_iterative[1][1]


def test_rpc_direct_iterative_nan():
    """
    test iterative direct localization with nan
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (col, row, alt) = (np.array([600, np.nan]), np.array([200, 210]), np.array([125]))
    p_direct0 = fctrat.direct_loc_inverse_iterative(row, col, alt, fill_nan=True)
    p_direct1 = fctrat.direct_loc_inverse_iterative(row[0], col[0], alt, fill_nan=True)
    #Point error col = 600, row = 200
    assert p_direct0[0][0] == p_direct1[0]
    assert p_direct0[1][0] == p_direct1[1]


def test_rpc_direct_iterative_all_nan():
    """
    test iterative direct localization with all nan
    """
    # comparasion between tabs en tab containing Nan cf. issue #46
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (col, row, alt) = (np.array([600, np.nan]), np.array([200, 210]), np.array([125]))
    direct_loc_tab = fctrat.direct_loc_inverse_iterative(row, col, alt, fill_nan=True)
    direct_loc_one_value = fctrat.direct_loc_inverse_iterative(row[0], col[0], alt, fill_nan=True)
    (col, row, alt) =(np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), np.array([125]))
    direct_loc_nan = fctrat.direct_loc_inverse_iterative(row, col, alt, 10, False)
    assert np.all(np.isnan(direct_loc_nan[0]))
    assert np.all(np.isnan(direct_loc_nan[1]))

    #Point error col = 600, row = 200
    assert direct_loc_tab[0][0] == direct_loc_one_value[0][0]
    assert direct_loc_tab[1][0] == direct_loc_one_value[1][0]
    assert direct_loc_tab[0][1] == fctrat.offset_X
    assert direct_loc_tab[1][1] == fctrat.offset_Y


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_direct_inverse_iterative(col, row, alt):
    """
    test the sequence of a direct localization followed by a inverse localization
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (lon, lat, __) = fctrat.direct_loc_h(row, col, alt)
    (row_inv, col_inv, __) = fctrat.inverse_loc(lon, lat, alt)
    (lon_iter, lat_iter, __) = fctrat.direct_loc_inverse_iterative(row_inv, col_inv, alt)
    assert lon == lon_iter
    assert lat == lat_iter


def test_rpc_minmax():
    """
    test get min and max altitude
    """
    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    fichier_dimap = os.path.join(data_folder,"rpc/PHRDIMAP_{}.XML".format(id_scene))
    fctrat = RPC.from_any(fichier_dimap)
    (h_min, h_max) = fctrat.get_alt_min_max()
    assert h_min == 532.5
    assert h_max == 617.5
