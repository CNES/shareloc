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
Test module for localisation class shareloc/bindings/dtm_intersection.cpp
"""
# pylint: disable=duplicate-code
# Standard imports
import os

import numpy as np

# Third party imports
import pytest

import bindings_cpp
from shareloc.dtm_reader import dtm_reader

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection

# Shareloc test imports
from ..helpers import data_path


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "roi,roi_phys,fill_nodata,fill_value",  # add read_data arg once #293 finished
    [
        (None, False, None, 0.0),
        ([256, 256, 512, 512], False, None, 0.0),
        ([44.57333333333334, 5.426666666666667, 44.36, 5.64], True, None, 0.0),
        (None, False, "mean", 12.0),
    ],
)
def test_constructor(roi, roi_phys, fill_nodata, fill_value):  # TODO : add more constructor use cases  before merge
    """
    Test DTMIntersection cpp constructor
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        read_data=True,
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
@pytest.mark.parametrize(  # Add pos_row=nb_rows
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
        read_data=True,
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
        read_data=True,
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
        read_data=True,
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
