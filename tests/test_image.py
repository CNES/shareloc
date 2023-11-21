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
Test module for image class shareloc/image.py
"""
# Standard imports
import os

# Third party imports
import pytest

# Shareloc imports
from shareloc.image import Image
from shareloc.proj_utils import transform_index_to_physical_point, transform_physical_point_to_index

# Shareloc test imports
from .helpers import data_path


# pylint: disable=too-many-arguments
@pytest.mark.parametrize(
    "image_name, row,col, origin_row, origin_col, pixel_size_row, pixel_size_col, "
    "pixel_rotation_row, pixel_rotation_col",
    [
        ("right_image.tif", 100, 200.5, 5162, 4915.0, 1.0, 1.0, 0.0, 0.0),
        ("right_image_resample.tif", 100.0, 200.0, 5162.0, 4915.0, 2.0, 0.5, 0.0, 0.0),
        ("left_image_shear.tif", 100.0, 200.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.36397023426620234),
    ],
)
@pytest.mark.unit_tests
def test_image_metadata(
    image_name, row, col, origin_row, origin_col, pixel_size_row, pixel_size_col, pixel_rotation_row, pixel_rotation_col
):
    """
    Test image class
    """
    data_folder = data_path()
    image_filename = os.path.join(data_folder, "image/phr_ventoux/", image_name)

    my_image = Image(image_filename)
    assert my_image.origin_row == origin_row
    assert my_image.origin_col == origin_col
    assert my_image.pixel_size_row == pixel_size_row
    assert my_image.pixel_size_col == pixel_size_col
    assert my_image.pixel_rotation_row == pixel_rotation_row
    assert my_image.pixel_rotation_col == pixel_rotation_col

    [phys_row, phys_col] = transform_index_to_physical_point(my_image.transform, row, col)
    assert phys_row == origin_row + (row + 0.5) * pixel_size_row + (col + 0.5) * pixel_rotation_row
    assert phys_col == origin_col + (col + 0.5) * pixel_size_col + (row + 0.5) * pixel_rotation_col

    row_index, col_index = transform_physical_point_to_index(my_image.trans_inv, phys_row, phys_col)
    assert row == row_index
    assert col == col_index


# pylint: disable=too-many-arguments
@pytest.mark.parametrize(
    "image_name, row,col, origin_row, origin_col, pixel_size_row, pixel_size_col, "
    "pixel_rotation_row, pixel_rotation_col",
    [
        ("right_image.tif", 100, 200.5, 5162, 4915.0, 1.0, 1.0, 0.0, 0.0),
        ("right_image_resample.tif", 100.0, 200.0, 5162.0, 4915.0, 2.0, 0.5, 0.0, 0.0),
        ("left_image_shear.tif", 100.0, 200.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.36397023426620234),
    ],
)
@pytest.mark.unit_tests
def test_image_metadata_roi(
    image_name, row, col, origin_row, origin_col, pixel_size_row, pixel_size_col, pixel_rotation_row, pixel_rotation_col
):
    """
    Test image class with roi
    """

    data_folder = data_path()
    image_filename = os.path.join(data_folder, "image/phr_ventoux/", image_name)

    # roi in index
    start_row = 10
    start_col = 20
    end_row = 310
    end_col = 320
    my_image_roi = Image(image_filename, roi=[start_row, start_col, end_row, end_col])
    assert my_image_roi.pixel_size_row == pixel_size_row
    assert my_image_roi.pixel_size_col == pixel_size_col
    assert my_image_roi.pixel_rotation_row == pixel_rotation_row
    assert my_image_roi.pixel_rotation_col == pixel_rotation_col

    phy_valid_row = origin_row + (row + 0.5) * pixel_size_row + (col + 0.5) * pixel_rotation_row
    phy_valid_col = origin_col + (col + 0.5) * pixel_size_col + (row + 0.5) * pixel_rotation_col
    [phys_row, phys_col] = transform_index_to_physical_point(my_image_roi.transform, row - start_row, col - start_col)
    assert phys_row == phy_valid_row
    assert phys_col == phy_valid_col

    row_index, col_index = transform_physical_point_to_index(my_image_roi.trans_inv, phys_row, phys_col)
    assert row == row_index + start_row
    assert col == col_index + start_col

    # roi in physical space
    start_phys_row = origin_row + start_row * pixel_size_row + start_col * pixel_rotation_row
    start_phys_col = origin_col + start_col * pixel_size_col + start_row * pixel_rotation_col
    end_phys_row = origin_row + end_row * pixel_size_row + end_col * pixel_rotation_row
    end_phys_col = origin_col + end_col * pixel_size_col + end_row * pixel_rotation_col
    roi = [start_phys_row, start_phys_col, end_phys_row, end_phys_col]
    my_image_roi_phys = Image(image_filename, roi=roi, roi_is_in_physical_space=True)
    [phys_row, phys_col] = transform_index_to_physical_point(
        my_image_roi_phys.transform, row - start_row, col - start_col
    )
    assert phys_row == phy_valid_row
    assert phys_col == phy_valid_col
