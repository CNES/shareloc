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
Test module for localisation class shareloc/localisation.py
"""


import os
import pytest
import numpy as np
from utils import test_path

from shareloc.dtm import DTM


@pytest.mark.parametrize("index_col,index_row, valid_alt", [(10.0, 20.0, 196.0), (20.5, 25.5, 189.5)])
@pytest.mark.unit_tests
def test_interp_dtm(index_col, index_row, valid_alt):
    """
    Test interp function
    """
    data_folder = test_path(alti="geoide", id_scene="P1BP--2017030824934340CP")
    # chargement du mnt
    fic = os.path.join(data_folder, "MNT_extrait/mnt_extrait.c1")
    dtmbsq = DTM(fic)
    vect_index = [index_row, index_col]
    coords = dtmbsq.index_to_ter(vect_index)
    lon_ref = 57.2083333333
    lat_ref = 22.2166666667
    pas_lon = 0.00833333333333
    pas_lat = -0.00833333333333
    assert dtmbsq.dtm_image.datum == "geoid"
    assert coords[0] == pytest.approx(lon_ref + index_col * pas_lon, 1e-12)
    assert coords[1] == pytest.approx(lat_ref + index_row * pas_lat, 1e-12)
    index = dtmbsq.ter_to_index(coords)
    assert index[1] == pytest.approx(index_col, 1e-12)
    assert index[0] == pytest.approx(index_row, 1e-12)
    alti = dtmbsq.interpolate(index_row, index_col)
    assert alti == valid_alt


#  geoid value is 50.73
@pytest.mark.parametrize("index_lon,index_lat, valid_alt", [(5.002777777777778, 44.99444444444445, 221.73)])
@pytest.mark.unit_tests
def test_interp_dtm_geoid(index_lon, index_lat, valid_alt):
    """
    Test interp function
    """
    dtm_file = os.path.join(os.environ["TESTPATH"], "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(os.environ["TESTPATH"], "dtm", "geoid", "egm96_15.gtx")
    dtm_ventoux = DTM(dtm_file, geoid_file)
    coords = [index_lon, index_lat]
    index = dtm_ventoux.ter_to_index(coords)
    alti = dtm_ventoux.interpolate(index[0], index[1])
    assert alti == pytest.approx(valid_alt, 1e-2)


@pytest.mark.unit_tests
def test_dtm_alt_min_max():
    """
    Test dtm alt min/max
    """
    dtm_file = os.path.join(os.environ["TESTPATH"], "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(os.environ["TESTPATH"], "dtm", "geoid", "egm96_15.gtx")
    dtm_ventoux = DTM(dtm_file, geoid_file)
    alt_min = dtm_ventoux.alt_min_cell
    alt_max = dtm_ventoux.alt_max_cell
    alt_valid_min = os.path.join(os.environ["TESTPATH"], "srtm_ventoux_alt_min.npy")
    alt_valid_max = os.path.join(os.environ["TESTPATH"], "srtm_ventoux_alt_max.npy")
    alt_min_vt = np.load(alt_valid_min)
    alt_max_vt = np.load(alt_valid_max)
    np.testing.assert_array_equal(alt_min, alt_min_vt)
    np.testing.assert_array_equal(alt_max, alt_max_vt)
