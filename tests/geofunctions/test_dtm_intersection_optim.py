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
Test module for localisation class shareloc/bindings/dtm_intersection.cpp
"""
# pylint: disable=duplicate-code
# Standard imports
import os
import pickle

import numpy as np

# Third party imports
import pytest

import bindings_cpp
from shareloc.dtm_reader import dtm_reader

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection
from shareloc.geomodels import GeoModel

# Shareloc test imports
from ..helpers import data_path


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "roi,roi_phys,fill_nodata,fill_value",
    [
        (None, False, None, 0.0),
        ([256, 256, 512, 512], False, None, 0.0),
        ([44.57333333333334, 5.426666666666667, 44.36, 5.64], True, None, 0.0),
        ([44.57333333333334, 5.426666666666667, 43.36, 6.94], True, None, 0.0),
        (None, False, "mean", 12.0),
    ],
)
def test_constructor(roi, roi_phys, fill_nodata, fill_value):
    """
    Test DTMIntersection cpp constructor
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=roi,
        roi_is_in_physical_space=roi_phys,
        fill_nodata=fill_nodata,
        fill_value=fill_value,
    )
    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    alt_min_cell_optim = np.array(dtm_ventoux_optim.get_alt_min_cell())
    alt_max_cell_optim = np.array(dtm_ventoux_optim.get_alt_max_cell())

    np.testing.assert_array_equal(dtm_ventoux_py.alt_min_cell.flatten(), alt_min_cell_optim)
    np.testing.assert_array_equal(dtm_ventoux_py.alt_max_cell.flatten(), alt_max_cell_optim)

    np.testing.assert_array_equal(dtm_ventoux_py.alt_data.flatten(), np.array(dtm_ventoux_optim.get_alt_data()))

    np.testing.assert_array_equal(dtm_ventoux_py.plane_coef_a.flatten(), np.array(dtm_ventoux_optim.get_plane_coef_a()))
    np.testing.assert_array_equal(dtm_ventoux_py.plane_coef_b.flatten(), np.array(dtm_ventoux_optim.get_plane_coef_b()))
    np.testing.assert_array_equal(dtm_ventoux_py.plane_coef_c.flatten(), np.array(dtm_ventoux_optim.get_plane_coef_c()))
    np.testing.assert_array_equal(dtm_ventoux_py.plane_coef_d.flatten(), np.array(dtm_ventoux_optim.get_plane_coef_d()))

    np.testing.assert_array_equal(dtm_ventoux_py.plans.flatten(), np.array(dtm_ventoux_optim.get_plans()))

    assert dtm_ventoux_py.alt_min == dtm_ventoux_optim.get_alt_min()
    assert dtm_ventoux_py.alt_max == dtm_ventoux_optim.get_alt_max()

    assert dtm_ventoux_py.tol_z == dtm_ventoux_optim.get_tol_z()
    assert dtm_ventoux_py.epsg == dtm_ventoux_optim.get_epsg()

    assert dtm_ventoux_py.nb_rows == dtm_ventoux_optim.get_nb_rows()
    assert dtm_ventoux_py.nb_columns == dtm_ventoux_optim.get_nb_columns()

    np.testing.assert_array_equal(
        np.array(dtm_ventoux_py.trans_inv.to_gdal()), np.array(dtm_ventoux_optim.get_trans_inv())
    )
    np.testing.assert_array_equal(
        np.array(dtm_ventoux_py.transform.to_gdal()), np.array(dtm_ventoux_optim.get_transform())
    )


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "i,positon",
    [
        (0, [5.00, 44.0, -100]),
        (1, [5.15, 44.15, 0]),
        (2, [5.35, 44.35, 50]),
        (3, [5.65, 44.65, 500]),
        (4, [5.85, 44.85, 200]),
        (5, [6.0, 45.0, 12000]),
    ],
)
def test_eq_plan(i, positon):
    """
    Test DTMIntersection eq_plan methode
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=[256, 256, 512, 512],
        roi_is_in_physical_space=False,
        fill_nodata=None,
        fill_value=0.0,
    )
    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    assert dtm_ventoux_py.eq_plan(i, positon) == dtm_ventoux_optim.eq_plan(i, positon)


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "pos_row,pos_col",
    [
        (0, 0),
        (1201, 0),
        (0, 1201),
        (450.236, 25.39),
        (2.369, 658.36),
        (1201, 1201),
    ],
)
def test_interpolate(pos_row, pos_col):
    """
    Test DTMIntersection interpolate methode
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=[256, 256, 512, 512],
        roi_is_in_physical_space=False,
        fill_nodata=None,
        fill_value=0.0,
    )
    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    res_py = dtm_ventoux_py.interpolate(pos_row, pos_col)
    res_optim = dtm_ventoux_optim.interpolate(pos_row, pos_col)

    assert res_py == res_optim


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "coords",
    [
        ([0, 0, 0]),
        ([1201, 0, 0]),
        ([0, 1201, 0]),
        ([450.236, 25.39, 0]),
        ([2.369, 658.36, 0]),
        ([1201, 1201, 0]),
    ],
)
def test_index_ter_methodes(coords):
    """
    Test DTMIntersection index_to_ter and ter_to_index methode
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=[256, 256, 512, 512],
        roi_is_in_physical_space=False,
        fill_nodata=None,
        fill_value=0.0,
    )
    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    res_py = dtm_ventoux_py.index_to_ter(coords)
    res_cpp = dtm_ventoux_optim.index_to_ter(coords)
    np.testing.assert_array_equal(res_py, np.array(res_cpp))

    coords_py = dtm_ventoux_py.ter_to_index(res_py)
    coords_cpp = dtm_ventoux_optim.ter_to_index(res_py)
    np.testing.assert_array_equal(coords_py, np.array(coords_cpp))

    # Identity
    np.testing.assert_allclose(np.array(coords), np.array(coords_cpp), 0, 7e-12)
    np.testing.assert_allclose(np.array(coords), np.array(coords_py), 0, 7e-12)


def test_intersect_dtm_cube():
    """
    test intersect_dtm_cube methode
    """

    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=None,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    )

    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    # -- No solution --

    # los from test_sensor_loc_dir_dtm -> not in adequation with the dtm
    los = np.array(
        [
            [57.37765525, 22.13204304, 4900.0],
            [57.37805121, 22.1312895, 3667.5],
            [57.37844733, 22.13053566, 2435.0],
            [57.37884359, 22.12978153, 1202.5],
            [57.37924, 22.12902711, -30.0],
        ]
    )

    (b_trouve_py, position_cube_py, alti_py, los_index_py) = dtm_ventoux_py.intersect_dtm_cube(los)

    los_x = los[:, 0].tolist()
    los_y = los[:, 1].tolist()
    los_z = los[:, 2].tolist()

    (b_trouve_cpp, position_cube_cpp, alti_cpp, los_index_x_cpp, los_index_y_cpp, los_index_z_cpp) = (
        dtm_ventoux_optim.intersect_dtm_cube(los_x, los_y, los_z)
    )

    assert b_trouve_py == b_trouve_cpp

    assert position_cube_py is None
    assert position_cube_cpp == [0.0, 0.0, 0.0]

    assert alti_py is None
    assert alti_cpp == 0.0

    np.testing.assert_array_equal(los_index_py[:, 0], los_index_x_cpp)
    np.testing.assert_array_equal(los_index_py[:, 1], los_index_y_cpp)
    np.testing.assert_array_equal(los_index_py[:, 2], los_index_z_cpp)

    # -- SOLUTION --

    los = np.array(
        [
            [
                5.19490591e00,
                4.42110907e01,
                2.81029295e03,
            ],
            [5.17180698e00, 4.41643824e01, -3.27181677e04],
        ]
    )

    (b_trouve_py, position_cube_py, alti_py, los_index_py) = dtm_ventoux_py.intersect_dtm_cube(los)

    los_x = los[:, 0].tolist()
    los_y = los[:, 1].tolist()
    los_z = los[:, 2].tolist()

    (b_trouve_cpp, position_cube_cpp, alti_cpp, los_index_x_cpp, los_index_y_cpp, los_index_z_cpp) = (
        dtm_ventoux_optim.intersect_dtm_cube(los_x, los_y, los_z)
    )

    assert b_trouve_py == b_trouve_cpp

    np.testing.assert_array_equal(position_cube_py, position_cube_cpp)

    assert alti_py == alti_cpp

    np.testing.assert_array_equal(los_index_py[:, 0], los_index_x_cpp)
    np.testing.assert_array_equal(los_index_py[:, 1], los_index_y_cpp)
    np.testing.assert_array_equal(los_index_py[:, 2], los_index_z_cpp)


def test_init_min_max():
    """
    test init_min_max fonction
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=None,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    )

    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    alt_data = np.array(dtm_ventoux_py.alt_data).flatten()
    res_cpp1, res_cpp2 = bindings_cpp.init_min_max(alt_data, dtm_image.nb_rows, dtm_image.nb_columns)

    np.testing.assert_array_equal(res_cpp1, dtm_ventoux_py.alt_min_cell.flatten())
    np.testing.assert_array_equal(res_cpp2, dtm_ventoux_py.alt_max_cell.flatten())


def test_intersection():
    """
    test intersection methode
    """

    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=None,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    )

    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    point_b = [946.6927376111748, 233.88631181467858, 2809.292947565206]
    h_intersect = 2.8146517369487256e-05
    los_x_index = [946.6911600000021, 1002.7411199999988]
    los_y_index = [233.88709199999994, 206.16837600000053]
    los_z_index = [2810.29295, -32718.1677]

    (bool_2_py, point_r_py) = dtm_ventoux_py.intersection(
        np.array([los_x_index, los_y_index, los_z_index]).T, np.array(point_b), h_intersect
    )

    (bool_2_cpp, point_r_x_cpp, point_r_y_cpp, point_r_z_cpp) = dtm_ventoux_optim.intersection(
        los_x_index, los_y_index, los_z_index, point_b, h_intersect
    )

    assert bool_2_py == bool_2_cpp
    assert point_r_py[0] == point_r_x_cpp
    assert point_r_py[1] == point_r_y_cpp
    assert point_r_py[2] == point_r_z_cpp


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "roi,roi_phys,fill_nodata,fill_value",
    [
        (None, False, None, 0.0),
        ([256, 256, 512, 512], False, None, 0.0),
        ([44.57333333333334, 5.426666666666667, 44.36, 5.64], True, None, 0.0),
        (None, False, "mean", 12.0),
    ],
)
def test_serialization(roi, roi_phys, fill_nodata, fill_value):
    """
    Test DTMIntersection serialization
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=roi,
        roi_is_in_physical_space=roi_phys,
        fill_nodata=fill_nodata,
        fill_value=fill_value,
    )

    dtm_original = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    data = pickle.dumps(dtm_original)
    dtm_deserialized = pickle.loads(data)

    np.testing.assert_array_equal(
        np.array(dtm_deserialized.get_alt_min_cell()), np.array(dtm_original.get_alt_min_cell())
    )
    np.testing.assert_array_equal(
        np.array(dtm_deserialized.get_alt_max_cell()), np.array(dtm_original.get_alt_max_cell())
    )

    np.testing.assert_array_equal(np.array(dtm_deserialized.get_alt_data()), np.array(dtm_original.get_alt_data()))

    np.testing.assert_array_equal(
        np.array(dtm_deserialized.get_plane_coef_a()), np.array(dtm_original.get_plane_coef_a())
    )
    np.testing.assert_array_equal(
        np.array(dtm_deserialized.get_plane_coef_b()), np.array(dtm_original.get_plane_coef_b())
    )
    np.testing.assert_array_equal(
        np.array(dtm_deserialized.get_plane_coef_c()), np.array(dtm_original.get_plane_coef_c())
    )
    np.testing.assert_array_equal(
        np.array(dtm_deserialized.get_plane_coef_d()), np.array(dtm_original.get_plane_coef_d())
    )

    np.testing.assert_array_equal(np.array(dtm_deserialized.get_plans()), np.array(dtm_original.get_plans()))

    assert dtm_deserialized.get_alt_min() == dtm_original.get_alt_min()
    assert dtm_deserialized.get_alt_max() == dtm_original.get_alt_max()

    assert dtm_deserialized.get_tol_z() == dtm_original.get_tol_z()
    assert dtm_deserialized.get_epsg() == dtm_original.get_epsg()

    assert dtm_deserialized.get_nb_rows() == dtm_original.get_nb_rows()
    assert dtm_deserialized.get_nb_columns() == dtm_original.get_nb_columns()

    np.testing.assert_array_equal(np.array(dtm_deserialized.get_trans_inv()), np.array(dtm_original.get_trans_inv()))
    np.testing.assert_array_equal(np.array(dtm_deserialized.get_transform()), np.array(dtm_original.get_transform()))

    # Test intersection

    rpc_file = os.path.join(data_path(), "rpc/phr_ventoux/RPC_PHR1B_P_201308051042194_SEN_690908101-001.XML")
    rpc_optim = GeoModel(rpc_file, "RPCoptim")

    rows = [i * 10 for i in range(10)]
    cols = [i * 10 for i in range(10)]

    res_original = rpc_optim.direct_loc_dtm(rows, cols, dtm_original)
    res_deserialized = rpc_optim.direct_loc_dtm(rows, cols, dtm_deserialized)

    np.testing.assert_array_equal(res_original, res_deserialized)


def test_intersection_n_los_dtm():
    """
    Test intersection_n_los_dtm method
    """
    mnt = os.path.join(data_path(), "dtm/srtm_ventoux/srtm90_non_void_filled/N44E005.hgt")
    dtm_image = dtm_reader(mnt)

    dtm_ventoux_py = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    dtm_ventoux_optim = bindings_cpp.DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    # los from rpc.py direct_loc_dtm of test_direct_loc_dtm from test_rpc_optim
    los = np.array(
        [
            [[5.16280369, 44.23317309, 2758.0], [5.16086044, 44.22958834, 31.0]],
            [[5.16284328, 44.23142179, 2758.0], [5.16090279, 44.22783679, 31.0]],
            [[5.16288289, 44.22967056, 2758.0], [5.16094517, 44.22608531, 31.0]],
        ]
    )

    res_optim = dtm_ventoux_optim.intersection_n_los_dtm(los)
    res_py = dtm_ventoux_py.intersection_n_los_dtm(los)

    np.testing.assert_array_equal(res_optim, res_py)
