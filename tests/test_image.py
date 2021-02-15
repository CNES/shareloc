#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2021 Centre National d'Etudes Spatiales (CNES).
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
Test module for image class shareloc/image/image.py
"""

import os
import pytest
from utils import test_path
from shareloc.image.image import Image


@pytest.mark.parametrize(
    "image_name, row,col, origin_row, origin_col, pixel_size_row, pixel_size_col",
    [
        ("right_image.tif", 100, 200.5, 5162, 4915.0, 1.0, 1.0),
        ("right_image_resample.tif", 100.0, 200.0, 5162.0, 4915.0, 2.0, 0.5),
    ],
)
@pytest.mark.unit_tests
def test_image_metadata(image_name, row, col, origin_row, origin_col, pixel_size_row, pixel_size_col):
    """
    Test image class
    """
    data_folder = test_path()
    image_filename = os.path.join(data_folder, "image/phr_ventoux/", image_name)

    my_image = Image(image_filename)
    assert my_image.origin_row == origin_row
    assert my_image.origin_col == origin_col
    assert my_image.pixel_size_row == pixel_size_row
    assert my_image.pixel_size_col == pixel_size_col

    [phys_row, phys_col] = my_image.transform_index_to_physical_point(row, col)
    assert phys_row == origin_row + (row + 0.5) * pixel_size_row
    assert phys_col == origin_col + (col + 0.5) * pixel_size_col

    row_index, col_index = my_image.transform_physical_point_to_index(phys_row, phys_col)
    assert row == row_index
    assert col == col_index
