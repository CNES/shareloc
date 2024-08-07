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
Test module for localisation class shareloc/geofunctions/dtm_intersection.py
"""
# pylint: disable=duplicate-code
# Standard imports
import os

import numpy as np

# Third party imports
import pytest

from shareloc.dtm_reader import dtm_reader, interpolate_geoid_height

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection

# Shareloc test imports
from ..helpers import data_path


@pytest.mark.parametrize(
    "lon,lat, valid_alt",
    [
        (0.1, 0.0, 17.12881622),
        (0.0, 0.0, 17.16157913),
        (10.125, 0.0, 8.72032166),
        (5.19368066, 44.20749145, 50.8618354),
        (5.17381589, 44.22559175, 50.87610),
        (-119.0, -9, -12.516327),
        (179.875, 44.0, -7.975677),
        (np.nan, np.nan, np.nan),
    ],
)
@pytest.mark.unit_tests
def test_geoid_height(lon, lat, valid_alt):
    """
    Test interpolate geoid height
    """
    geoid_file = os.path.join(data_path(), "dtm/geoid/egm96_15.gtx")
    positions = np.asarray([lon, lat])[np.newaxis, :]
    geoid_height = interpolate_geoid_height(geoid_file, positions)[0]
    assert geoid_height == pytest.approx(valid_alt, abs=1e-6, nan_ok=True)

    geoid_file = os.path.join(data_path(), "dtm/geoid/egm96.grd")
    geoid_height = interpolate_geoid_height(geoid_file, positions)[0]
    # ground truth have been computed using .gtx small alt diff between .gtx and .grd
    assert geoid_height == pytest.approx(valid_alt, abs=1e-3, nan_ok=True)


@pytest.mark.parametrize("index_col,index_row, valid_alt", [(10.0, 20.0, 196.0), (20.5, 25.5, 189.5)])
@pytest.mark.unit_tests
def test_interp_dtm(index_col, index_row, valid_alt):
    """
    Test interp function
    """
    data_folder = data_path("geoide", "P1BP--2017030824934340CP")
    # load MNT
    fic = os.path.join(data_folder, "MNT_P1BP--2017030824934340CP.tif")
    dtm_image = dtm_reader(
        fic,
        None,
        roi=None,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    )
    dtm = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    vect_index = [index_row, index_col]
    coords = dtm.index_to_ter(vect_index)
    lon_ref = 57.2083333333
    lat_ref = 22.2166666667
    pas_lon = 0.00833333333333
    pas_lat = -0.00833333333333
    assert coords[0] == pytest.approx(lon_ref + index_col * pas_lon, 1e-12)
    assert coords[1] == pytest.approx(lat_ref + index_row * pas_lat, 1e-12)
    index = dtm.ter_to_index(coords)
    assert index[1] == pytest.approx(index_col, 1e-12)
    assert index[0] == pytest.approx(index_row, 1e-12)
    alti = dtm.interpolate(index_row, index_col)
    assert alti == valid_alt


#  geoid value is 50.73
@pytest.mark.parametrize("index_lon,index_lat, valid_alt", [(5.002777777777778, 44.99444444444445, 221.73)])
@pytest.mark.unit_tests
def test_interp_dtm_geoid(index_lon, index_lat, valid_alt):
    """
    Test interp function
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
    dtm_ventoux = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    coords = [index_lon, index_lat]
    index = dtm_ventoux.ter_to_index(coords)
    alti = dtm_ventoux.interpolate(index[0], index[1])
    assert alti == pytest.approx(valid_alt, 1e-2)


@pytest.mark.unit_tests
def test_dtm_alt_min_max():
    """
    Test dtm alt min/max
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
    dtm_ventoux = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    alt_min = dtm_ventoux.alt_min_cell
    alt_max = dtm_ventoux.alt_max_cell
    alt_valid_min = os.path.join(data_path(), "srtm_ventoux_alt_min.npy")
    alt_valid_max = os.path.join(data_path(), "srtm_ventoux_alt_max.npy")
    alt_min_vt = np.load(alt_valid_min)
    alt_max_vt = np.load(alt_valid_max)
    np.testing.assert_array_equal(alt_min, alt_min_vt)
    np.testing.assert_array_equal(alt_max, alt_max_vt)


@pytest.mark.unit_tests
@pytest.mark.parametrize(
    "roi,data_size_ref",
    [
        ([-5, -5, 256, 256], (256, 256)),
        ([0, 0, 1201, 1201], (1201, 1201)),  # dtm size
        ([-15, -15, 1215, 1215], (1201, 1201)),
        ([5, 5, 100, 100], (95, 95)),
        ([5, 5, 1215, 1215], (1196, 1196)),
    ],
)
def test_init_dtm_reader_roi(roi, data_size_ref):
    """
    Test dtm reader with roi larger than image
    """
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=roi,
        roi_is_in_physical_space=False,
        fill_nodata=None,
        fill_value=0.0,
    )
    _ = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    assert dtm_image.alt_data.shape == data_size_ref

    # ROI in physical space
    footprint_dtm = [43.9995833333333337, 4.9995833333333337, 45.0004166666666663, 6.0004166666666672]
    footprint_dtm[0] -= 1
    footprint_dtm[1] -= 1
    footprint_dtm[2] += 1
    footprint_dtm[3] += 1

    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=footprint_dtm,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    )
    _ = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    assert dtm_image.alt_data.shape == (1201, 1201)
