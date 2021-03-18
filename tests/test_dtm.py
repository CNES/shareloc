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
from utils import test_path

from shareloc.dtm import DTM


@pytest.mark.parametrize("index_x,index_y", [(10.5, 20.5)])
@pytest.mark.parametrize("valid_alt", [(198.0)])
@pytest.mark.unit_tests
def test_interp_dtm(index_x, index_y, valid_alt):
    """
    Test interp function
    """
    data_folder = test_path(alti="geoide", id_scene="P1BP--2017030824934340CP")
    # chargement du mnt
    fic = os.path.join(data_folder, "MNT_extrait/mnt_extrait.c1")
    dtmbsq = DTM(fic)
    vect_index = [index_x, index_y]
    coords = dtmbsq.index_to_ter(vect_index)
    lon_ref = 57.2083333333
    lat_ref = 22.2166666667
    pas_lon = 0.00833333333333
    pas_lat = -0.00833333333333
    assert coords[0] == lon_ref + index_y * pas_lon
    assert coords[1] == lat_ref + index_x * pas_lat
    index = dtmbsq.ter_to_index(coords)
    assert index[0] == pytest.approx(index_x, 1e-12)
    assert index[1] == pytest.approx(index_y, 1e-12)
    alti = dtmbsq.interpolate(index_x - 0.5, index_y - 0.5)
    assert alti == valid_alt
