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
import os
import pytest
from utils import test_path

from shareloc.localization import Localization
from shareloc.image.image  import Image

@pytest.mark.parametrize("row,col", [(100.0, 200.0)])
@pytest.mark.unit_tests
def test_image_metadata(row, col):
    """
     Test image class
    """
    data_folder = test_path()
    image_filename = os.path.join(data_folder, 'image/phr_ventoux/right_image.tif')

    my_image = Image(image_filename)
    assert my_image.origin_row == 5162.0
    assert my_image.origin_col == 4915.0
    assert my_image.pixel_size_row == 1
    assert my_image.pixel_size_col == 1

    [phys_row, phys_col] = my_image.transform_index_to_physical_point(row, col)

    assert phys_row == my_image.origin_row + (row + 0.5) * my_image.pixel_size_row
    assert phys_col == my_image.origin_col + (col + 0.5) * my_image.pixel_size_col