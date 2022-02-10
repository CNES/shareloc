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
Test module for proj utils  function  shareloc/proj_utils.py
"""

# Third party imports
import numpy as np
import pytest

# Shareloc imports
from shareloc.proj_utils import coordinates_conversion


@pytest.mark.unit_tests
def test_coordinates_conversion():
    """
    Test coordinates conversion
    """
    in_crs = 4326
    out_crs = 4978

    point_wgs84 = np.asarray([[7.05396752, 43.73000865, 4900.0], [7.05860411, 43.72347311, -30.0]])
    point_ecef = coordinates_conversion(point_wgs84, in_crs, out_crs)
    coords_vt_ecef = np.asarray(
        [[4584837.334948, 567331.361674, 4389850.562378], [4581754.08394, 567326.291517, 4385917.904472]]
    )
    np.testing.assert_allclose(point_ecef, coords_vt_ecef, atol=1e-5, rtol=0)
