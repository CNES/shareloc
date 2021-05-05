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
Test module for euclidium utils functions
"""


import pytest


from shareloc.euclidium_utils import identify_gdlib_code


@pytest.mark.unit_tests
def test_identify_gdlib_code():
    """
    Test identify gdlib code
    """
    srs, datum = identify_gdlib_code("WGS84:G-D/WGS84:Z-M")
    assert srs == 4326
    assert datum == "ellipsoid"
    srs, datum = identify_gdlib_code("GRS80:G-D/:H-M")
    assert srs == 4269
    assert datum == "geoid"
    srs, datum = identify_gdlib_code("32631/:Z-M")
    assert srs == 32631
    assert datum == "ellipsoid"
