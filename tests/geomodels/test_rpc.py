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
Module to test functions that use rpc
"""
# TODO: refactor to disable no-member in RPC class
# pylint: disable=no-member


# Standard imports
import os

# Third party imports
import numpy as np
import pytest
import rasterio

from shareloc.dtm_reader import dtm_reader

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection
from shareloc.geomodels import GeoModel
from shareloc.geomodels.rpc_readers import identify_dimap, identify_ossim_kwl, rpc_reader_via_rasterio

# Shareloc test imports
from ..helpers import data_path


def test_rpc_drivers():
    """
    test rpc driver identification
    """
    data_folder = data_path()

    # Test DIMAP RPC
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")
    fctrat_dimap = GeoModel(file_dimap)
    assert fctrat_dimap.driver_type == "dimap_v1.4"

    # Test OSSIM KWL GEOM RPC
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    fctrat_geom = GeoModel(file_geom)
    assert fctrat_geom.driver_type == "ossim_kwl"

    # Test DIMAPv2 RPC
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_dimap_v2 = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    fctrat_dimap_v2 = GeoModel(file_dimap_v2)
    assert fctrat_dimap_v2.driver_type == "dimap_v2.15"

    # Test fake RPC
    fake_rpc = os.path.join(data_folder, "rpc/fake_rpc.txt")
    try:
        GeoModel(fake_rpc)  # Raise ValueError->True
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


@pytest.mark.parametrize(
    "id_scene,lon,lat,alt, col_vt,row_vt",
    [
        (
            "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001",
            7.048662660737769592,
            43.72774839443545858,
            0.0,
            100.5,
            200.5,
        ),
        (
            "PHR1B_P_201709281038045_SEN_PRG_FC_178608-001",
            7.11526088296757386331137240632,
            43.684281179313565246502548689,
            850.0,
            10121.5657,
            10236.4310,
        ),
    ],
)
def test_rpc_ossim_kwl(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization from ossim file
    """
    data_folder = data_path()
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    fctrat_geom = GeoModel(file_geom)

    (row, col, __) = fctrat_geom.inverse_loc(lon, lat, alt)
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)


@pytest.mark.parametrize(
    "id_scene,lon,lat,alt, col_vt,row_vt",
    [
        (
            "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.geom",
            7.048662660737769592,
            43.72774839443545858,
            0.0,
            100.5,
            200.5,
        ),
        (
            "RPC_PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.XML",
            7.048662660737769592,
            43.72774839443545858,
            0.0,
            100.5,
            200.5,
        ),
        (
            "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif",
            7.048662660737769592,
            43.72774839443545858,
            0.0,
            100.5,
            200.5,
        ),
    ],
)
def test_rpc_from_any_file(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization from any file
    """
    data_folder = data_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)
    fctrat = GeoModel(rpc_file)
    (row, col, __) = fctrat.inverse_loc(lon, lat, alt)
    assert fctrat.epsg == 4326
    assert fctrat.datum == "ellipsoid"
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)


@pytest.mark.parametrize(
    "id_scene,lon,lat,alt",
    [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif", 7.048662660737769592, 43.72774839443545858, 0.0)],
)
def test_rpc_direct_iterative(id_scene, lon, lat, alt):
    """
    test the sequence of a inverse localization followed by a direct localization using Geotif driver
    """
    data_folder = data_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)
    fctrat = GeoModel(rpc_file)
    (row, col, __) = fctrat.inverse_loc(lon, lat, alt)
    (lon2, lat2, __) = fctrat.direct_loc_inverse_iterative(row, col, alt)
    assert lon == pytest.approx(lon2, abs=1e-2)
    assert lat == pytest.approx(lat2, abs=1e-2)


@pytest.mark.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)
@pytest.mark.parametrize(
    "prod, can_read",
    [
        ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.tif", True),
        ("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001_nogeo.tif", False),
    ],
)
def test_rpc_from_geotiff_without_rpc(prod, can_read):
    """
    test file without rpc
    """
    data_folder = data_path()
    rpc_file = os.path.join(data_folder, "rpc", prod)
    try:
        GeoModel(rpc_file)
        assert can_read
    except ValueError:
        assert not can_read


@pytest.mark.parametrize(
    "id_scene,lon,lat,alt, col_vt,row_vt",
    [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001", 7.048662660737769592, 43.72774839443545858, 0.0, 100.5, 200.5)],
)
def test_rpc_dimap_v2(id_scene, lon, lat, alt, row_vt, col_vt):
    """
    test inverse localization using ossim and dimap files
    """
    data_folder = data_path()
    file_dimap = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    fctrat_dimap = GeoModel(file_dimap)

    (row, col, __) = fctrat_dimap.inverse_loc(lon, lat, alt)
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)

    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    fctrat_geom = GeoModel(file_geom)
    (row, col, __) = fctrat_geom.inverse_loc(lon, lat, alt)
    assert col == pytest.approx(col_vt, abs=1e-2)
    assert row == pytest.approx(row_vt, abs=1e-2)


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_phrdimap(col, row, alt):
    """
    test the sequence of a inverse localization followed by a direct localization using dimap file
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    fctrat = GeoModel(file_dimap)

    (lonlatalt) = fctrat.direct_loc_h(row, col, alt)

    (row_ar, col_ar, __) = fctrat.inverse_loc(lonlatalt[0][0], lonlatalt[0][1], lonlatalt[0][2])
    assert col_ar == pytest.approx(col, abs=2e-2)
    assert row_ar == pytest.approx(row, abs=2e-2)


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_direct_inverse_iterative_vs_direct(col, row, alt):
    """
    compare direct localization and iterative direct localization
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    fctrat = GeoModel(file_dimap)

    # (col,lig,alt)=(100,1000,400)
    lonlatalt = fctrat.direct_loc_h(row, col, alt)
    (x_inv, y_inv, __) = fctrat.direct_loc_inverse_iterative(row, col, alt)

    #  Comparison direct loc RPC / inverse loc RPC iterative
    # The errors observed are due to the differences between loc dir and loc inv RPC
    assert lonlatalt[0][0] == pytest.approx(x_inv, abs=10.0 / 111111000)
    assert lonlatalt[0][1] == pytest.approx(y_inv, abs=10.0 / 111111000)


def test_rpc_direct_inverse_iterative_vs_direct_multiple_points():
    """
    compare direct localization and iterative direct localization with multiple points
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    fctrat = GeoModel(file_dimap)

    (col, row, alt) = (np.array([600, 610]), np.array([200, 210]), np.array([125]))
    p_direct = fctrat.direct_loc_h(row, col, alt)
    (rowinv, colinv, __) = fctrat.inverse_loc(p_direct[:, 0], p_direct[:, 1], alt)
    p_direct_iterative = fctrat.direct_loc_inverse_iterative(rowinv, colinv, alt)

    # Comparison direct loc RPC / inverse loc RPC iterative
    # The errors observed are due to the differences between loc dir and loc inv RPC

    # Point error col = 600, row = 200
    assert p_direct[0, 0] == p_direct_iterative[0][0]
    assert p_direct[0, 1] == p_direct_iterative[1][0]

    # Point error col = 601, row = 201
    assert p_direct[1, 0] == p_direct_iterative[0][1]
    assert p_direct[1, 1] == p_direct_iterative[1][1]


def test_rpc_direct_iterative_nan():
    """
    test iterative direct localization with nan
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    fctrat = GeoModel(file_dimap)

    (col, row, alt) = (np.array([600, np.nan]), np.array([200, 210]), np.array([125]))
    p_direct0 = fctrat.direct_loc_inverse_iterative(row, col, alt, fill_nan=True)
    p_direct1 = fctrat.direct_loc_inverse_iterative(row[0], col[0], alt, fill_nan=True)
    # Point error col = 600, row = 200
    assert p_direct0[0][0] == p_direct1[0]
    assert p_direct0[1][0] == p_direct1[1]

    p_direct0 = fctrat.direct_loc_h(row, col, alt, fill_nan=True)
    p_direct1 = fctrat.direct_loc_h(row[0], col[0], alt, fill_nan=True)

    np.testing.assert_array_equal(p_direct0[0, :], p_direct1[0, :])


def test_rpc_direct_iterative_all_nan():
    """
    test iterative direct localization with all nan
    """
    # comparison between tabs en tab containing Nan cf. issue #46
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    fctrat = GeoModel(file_dimap)

    (col, row, alt) = (np.array([600, np.nan]), np.array([200, 210]), np.array([125]))
    direct_loc_tab = fctrat.direct_loc_inverse_iterative(row, col, alt, fill_nan=True)
    direct_loc_one_value = fctrat.direct_loc_inverse_iterative(row[0], col[0], alt, fill_nan=True)
    (col, row, alt) = (np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), np.array([125]))
    direct_loc_nan = fctrat.direct_loc_inverse_iterative(row, col, alt, 10, False)
    assert np.all(np.isnan(direct_loc_nan[0]))
    assert np.all(np.isnan(direct_loc_nan[1]))
    direct_loc_nan = fctrat.direct_loc_h(row, col, alt, False)
    assert np.all(np.isnan(direct_loc_nan[:, 0]))
    assert np.all(np.isnan(direct_loc_nan[:, 1]))

    # Point error col = 600, row = 200
    assert direct_loc_tab[0][0] == direct_loc_one_value[0][0]
    assert direct_loc_tab[1][0] == direct_loc_one_value[1][0]
    assert direct_loc_tab[0][1] == fctrat.offset_x
    assert direct_loc_tab[1][1] == fctrat.offset_y


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_rpc_direct_inverse_iterative(col, row, alt):
    """
    test the sequence of a direct localization followed by a inverse localization
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")

    fctrat = GeoModel(file_dimap)

    lonlatalt = fctrat.direct_loc_h(row, col, alt)
    (row_inv, col_inv, __) = fctrat.inverse_loc(lonlatalt[0][0], lonlatalt[0][1], lonlatalt[0][2])
    (lon_iter, lat_iter, __) = fctrat.direct_loc_inverse_iterative(row_inv, col_inv, alt)
    assert lonlatalt[0][0] == lon_iter
    assert lonlatalt[0][1] == lat_iter


@pytest.mark.parametrize(
    "id_scene, index_x,index_y", [("RPC_PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.XML", 10.5, 20.5)]
)
@pytest.mark.unit_tests
def test_rpc_direct_dtm(id_scene, index_x, index_y):
    """
    Test direct localization on DTMIntersection
    """
    vect_index = [index_x, index_y]

    data_folder = data_path()
    rpc_file = os.path.join(data_folder, "rpc", id_scene)
    fctrat = GeoModel(rpc_file)
    id_scene = "P1BP--2017092838284574CP"
    data_folder_mnt = data_path("ellipsoide", id_scene)

    fic = os.path.join(data_folder_mnt, f"MNT_{id_scene}.tif")
    dtm_image = dtm_reader(fic)
    dtm = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    [lon, lat] = dtm.index_to_ter(vect_index)
    alt = dtm.interpolate(index_x, index_y)

    row, col, alt = fctrat.inverse_loc(lon, lat, alt)
    lonlath = fctrat.direct_loc_dtm(row[0], col[0], dtm)
    assert lon == pytest.approx(lonlath[0][0], abs=1e-7)
    assert lat == pytest.approx(lonlath[0][1], abs=1e-7)
    assert alt == pytest.approx(lonlath[0][2], abs=1e-4)


@pytest.mark.parametrize(
    "id_scene, col,row",
    [("PHR1B_P_201709281038393_SEN_PRG_FC_178609-001", 100.5, 200.5)],
)
def test_rpc_los_extrapolation(id_scene, row, col):
    """
    test los extrapolation
    """
    data_folder = data_path()
    file_dimap = os.path.join(data_folder, f"rpc/RPC_{id_scene}.XML")
    fctrat = GeoModel(file_dimap)
    los_edges = fctrat.los_extrema(row, col)
    altmin = -10
    altmax = 2000
    extrapoled_los_edges = fctrat.los_extrema(row, col, altmin, altmax)
    diff = los_edges[0, :] - los_edges[1, :]
    lon_extrapoled_max = los_edges[1, 0] + diff[0] * (altmax - los_edges[1, 2]) / diff[2]
    assert lon_extrapoled_max == extrapoled_los_edges[0, 0]


def test_rpc_minmax():
    """
    test get min and max altitude
    """
    data_folder = data_path()
    id_scene = "P1BP--2018122638935449CP"
    fichier_dimap = os.path.join(data_folder, f"rpc/PHRDIMAP_{id_scene}.XML")
    fctrat = GeoModel(fichier_dimap)
    (h_min, h_max) = fctrat.get_alt_min_max()
    assert h_min == 532.5
    assert h_max == 617.5


def test_read_via_rasterio():
    """
    Test rpc read via rasterio
    """

    # -- Assert RPC constructor
    data_folder = data_path()
    rpc_path = os.path.join(data_folder, "rpc/wv3_20.NTF")

    rpc = GeoModel(rpc_path)

    assert rpc.driver_type == "rasterio_rpc"

    # -- Assert dict values
    rpc_params_ref = {}
    rpc_params_ref["offset_alt"] = 31
    rpc_params_ref["scale_alt"] = 501
    rpc_params_ref["offset_y"] = -34.5043
    rpc_params_ref["scale_y"] = 0.0531
    rpc_params_ref["den_row"] = [
        1,
        6.918796e-05,
        0.0004416544,
        1.091165e-05,
        -4.703357e-07,
        1.180023e-06,
        -2.818314e-07,
        4.24824e-05,
        -0.0001854623,
        5.946361e-05,
        -1.263282e-08,
        2.805569e-08,
        1.145311e-05,
        0,
        2.289254e-06,
        -0.001310256,
        1.146713e-07,
        0,
        4.505145e-06,
        0,
    ]
    rpc_params_ref["num_row"] = [
        0.002401507,
        -0.002429637,
        1.002863,
        -0.0007633858,
        -6.513516e-05,
        9.272476e-08,
        -1.018064e-05,
        -0.000491502,
        -0.002590553,
        1.890887e-07,
        1.184859e-06,
        -5.752598e-08,
        1.426182e-06,
        -1.451179e-07,
        4.540262e-05,
        0.0002559039,
        5.966063e-05,
        -2.931054e-08,
        1.022565e-06,
        -4.526981e-08,
    ]
    rpc_params_ref["offset_row"] = 17495
    rpc_params_ref["scale_row"] = 17996
    rpc_params_ref["offset_x"] = -58.6024
    rpc_params_ref["scale_x"] = 0.0803
    rpc_params_ref["den_col"] = [
        1,
        0.0007669113,
        0.0003657245,
        -0.0004263618,
        -8.254525e-06,
        -1.611481e-06,
        5.995759e-07,
        4.906291e-06,
        1.828526e-06,
        -1.831867e-06,
        -1.851723e-08,
        0,
        1.280374e-07,
        0,
        -1.671405e-07,
        3.795708e-07,
        0,
        -1.725435e-08,
        0,
        0,
    ]
    rpc_params_ref["num_col"] = [
        0.005014126,
        -1.021695,
        -0.0002265161,
        0.02071451,
        0.0003603364,
        -0.0003155448,
        0.0001644524,
        -0.004269895,
        -0.0004006558,
        6.408265e-06,
        1.348922e-06,
        -1.285029e-05,
        2.041723e-06,
        1.659585e-06,
        -1.835126e-05,
        -6.090936e-05,
        4.43883e-08,
        -2.871058e-06,
        -1.479638e-07,
        -3.507909e-08,
    ]
    rpc_params_ref["offset_col"] = 20749
    rpc_params_ref["scale_col"] = 21250
    rpc_params_ref["num_x"] = None
    rpc_params_ref["den_x"] = None
    rpc_params_ref["num_y"] = None
    rpc_params_ref["den_y"] = None
    rpc_params_ref["driver_type"] = "rasterio_rpc"

    rpc_params = rpc_reader_via_rasterio(rpc_path, topleftconvention=False)

    assert rpc_params_ref == rpc_params

    # -- Identity
    row = 501
    col = 501
    alt = 10
    lonlatalt = rpc.direct_loc_h(row, col, alt)
    (row_calc, col_calc, alt) = rpc.inverse_loc(lonlatalt[:, 0], lonlatalt[:, 1], lonlatalt[:, 2])

    assert row_calc[0] == pytest.approx(row, abs=1e-9)
    assert col_calc[0] == pytest.approx(col, abs=1e-9)
