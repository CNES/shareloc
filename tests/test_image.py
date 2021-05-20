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
import numpy as np
from utils import test_path
from shareloc.image.image import Image
from shareloc.image.dtm_image import DTMImage


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

    start_row = 10
    start_col = 20
    end_row = 310
    end_col = 320
    my_image_roi = Image(image_filename, roi=[start_row, start_col, end_row, end_col])
    assert my_image_roi.pixel_size_row == pixel_size_row
    assert my_image_roi.pixel_size_col == pixel_size_col

    [phys_row, phys_col] = my_image_roi.transform_index_to_physical_point(row - start_row, col - start_col)
    assert phys_row == origin_row + (row + 0.5) * pixel_size_row
    assert phys_col == origin_col + (col + 0.5) * pixel_size_col

    row_index, col_index = my_image_roi.transform_physical_point_to_index(phys_row, phys_col)
    assert row == row_index + start_row
    assert col == col_index + start_col

    ## roi in physical space
    start_phys_row = origin_row + start_row * pixel_size_row
    start_phys_col = origin_col + start_col * pixel_size_col
    end_phys_row = origin_row + end_row * pixel_size_row
    end_phys_col = origin_col + end_col * pixel_size_col
    roi = [start_phys_row, start_phys_col, end_phys_row, end_phys_col]
    my_image_roi_phys = Image(image_filename, roi=roi, roi_is_in_physical_space=True)

    [phys_row, phys_col] = my_image_roi_phys.transform_index_to_physical_point(row - start_row, col - start_col)
    assert phys_row == origin_row + (row + 0.5) * pixel_size_row
    assert phys_col == origin_col + (col + 0.5) * pixel_size_col


@pytest.mark.parametrize(
    "row,col, origin_row, origin_col, pixel_size_row, pixel_size_col",
    [(100, 200.5, 43.74583333336666, 7.012500000003335, -0.00833333333333, 0.00833333333333)],
)
@pytest.mark.unit_tests
def test_bsqimage_metadata(row, col, origin_row, origin_col, pixel_size_row, pixel_size_col):
    """
    Test  dtmimage class
    """
    data_folder = test_path("ellipsoide", "P1BP--2017092838319324CP")
    # chargement du mnt
    image_filename = os.path.join(data_folder, "MNT_extrait/mnt_extrait.c1")

    my_image = DTMImage(image_filename)
    assert my_image.origin_row == origin_row
    assert my_image.origin_col == origin_col
    assert my_image.pixel_size_row == pixel_size_row
    assert my_image.pixel_size_col == pixel_size_col

    [phys_row, phys_col] = my_image.transform_index_to_physical_point(row, col)
    assert phys_row == origin_row + (row + 0.5) * pixel_size_row
    assert phys_col == origin_col + (col + 0.5) * pixel_size_col

    row_index, col_index = my_image.transform_physical_point_to_index(phys_row, phys_col)
    assert row == pytest.approx(row_index, abs=1e-12)
    assert col == pytest.approx(col_index, abs=1e-12)
    assert my_image.datum == "ellipsoid"


def test_dtm_fillnodata():
    """
    Test dtm image fillnodata
    """
    dtm_file = os.path.join(os.environ["TESTPATH"], "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    my_image_with_nodata = DTMImage(dtm_file, read_data=True, fill_nodata=None)
    nodata_index = np.argwhere(my_image_with_nodata.mask == 0)[0]
    assert my_image_with_nodata.data[0, nodata_index[1], nodata_index[2]] == -32768
    my_image_rio_fillnodata = DTMImage(dtm_file, read_data=True, fill_nodata="rio_fillnodata")
    assert my_image_rio_fillnodata.data[0, nodata_index[1], nodata_index[2]] == 783
    my_image_mean = DTMImage(dtm_file, read_data=True, fill_nodata="mean")
    assert my_image_mean.data[0, nodata_index[1], nodata_index[2]] == 872
