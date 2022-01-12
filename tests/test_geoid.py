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
Test module for geoid function class shareloc/localisation.py
"""


import os
import pytest
import numpy as np
from helpers import data_path

from shareloc.geoid import interpolate_geoid_height


@pytest.mark.parametrize(
    "lon,lat, valid_alt",
    [
        (0.0, 0.0, 17.16157913),
        (10.125, 0.0, 8.72032166),
        (5.19368066, 44.20749145, 50.8618354),
        (5.17381589, 44.22559175, 50.87610),
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
    assert geoid_height == pytest.approx(valid_alt, abs=1e-6)
