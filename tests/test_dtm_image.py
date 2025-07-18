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
Test module for image class shareloc/dtm_image.py
"""

# Standard imports
import os

# Third party imports
import numpy as np
import pytest

# Shareloc imports
from shareloc.dtm_reader import dtm_reader

# Shareloc test imports
from .helpers import data_path


@pytest.mark.unit_tests
def test_dtm_fillnodata():
    """
    Test dtm image fillnodata
    """
    # Set test data SRTM dtm path
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")

    my_image_with_nodata = dtm_reader(dtm_file, fill_nodata=None)
    nodata_index = np.argwhere(my_image_with_nodata.mask == 0)[0]
    assert my_image_with_nodata.data[nodata_index[0], nodata_index[1]] == -32768

    my_image_rio_fillnodata = dtm_reader(dtm_file, fill_nodata="rio_fillnodata")
    assert my_image_rio_fillnodata.data[nodata_index[0], nodata_index[1]] == 783

    my_image_mean = dtm_reader(dtm_file, fill_nodata="mean")
    assert my_image_mean.data[nodata_index[0], nodata_index[1]] == 872

    my_image_min = dtm_reader(dtm_file, fill_nodata="min")
    assert my_image_min.data[nodata_index[0], nodata_index[1]] == 32

    my_image_max = dtm_reader(dtm_file, fill_nodata="max")
    assert my_image_max.data[nodata_index[0], nodata_index[1]] == 2757

    my_image_median = dtm_reader(dtm_file, fill_nodata="median")
    assert my_image_median.data[nodata_index[0], nodata_index[1]] == 840

    my_image_constant = dtm_reader(dtm_file, fill_nodata="constant", fill_value=100.0)
    assert my_image_constant.data[nodata_index[0], nodata_index[1]] == 100.0

    # Set test data with hole in srtm
    dtm_file_srtm_hole = os.path.join(data_path(), "dtm", "srtm_ventoux", "N44E005_big_hole.tif")
    my_image_fill_hole = dtm_reader(dtm_file_srtm_hole, fill_nodata="rio_fillnodata")
    assert my_image_fill_hole.data[403, 1119] == 32


def test_dtm_with_nan():
    """
    Test dtm image containing nan
    """

    # image containing one NaN and nodata value
    dtm_file = os.path.join(data_path(), "dtm", "dtm_philondenx_32630_nan.tif")

    dtm_with_nan = dtm_reader(dtm_file, fill_nodata="mean")
    assert dtm_with_nan.stats["min"] == pytest.approx(210.45277, abs=1e-5)
    assert dtm_with_nan.stats["max"] == pytest.approx(244.91403, abs=1e-5)
    assert dtm_with_nan.stats["mean"] == pytest.approx(221.15659, abs=1e-5)
    assert dtm_with_nan.stats["median"] == pytest.approx(221.11378, abs=1e-5)
    assert np.sum(np.isnan(dtm_with_nan.data[:, :])) == 0
    # image containing one NaN and no nodata value
    dtm_file = os.path.join(data_path(), "dtm", "dtm_philondenx_32630_nan_no_nodata.tif")
    dtm_with_nan = dtm_reader(dtm_file, fill_nodata="mean")
    assert dtm_with_nan.stats["min"] == 0.0
    assert dtm_with_nan.stats["max"] == pytest.approx(244.91403, abs=1e-5)
    assert dtm_with_nan.stats["mean"] == pytest.approx(179.98820, abs=1e-5)
    assert dtm_with_nan.stats["median"] == pytest.approx(220.37231, abs=1e-5)
    assert np.sum(np.isnan(dtm_with_nan.data[:, :])) == 0


def test_dtm_outside_roi():
    """
    Raise error when the dtm is outside the roi
    """

    # Set test data SRTM dtm path
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")

    roi = [31.099861111111114, 29.950138888888887, 31.149861111111115, 30.000138888888888]

    with pytest.raises(RuntimeError) as excinfo:
        _ = dtm_reader(dtm_file, roi=roi, roi_is_in_physical_space=True, fill_nodata=None)

    assert str(excinfo.value) == (
        "the roi bounds are "
        "[29.950138888888887, 31.099861111111114, 30.000138888888888, "
        "31.149861111111115] while the dtm bounds are [4.999583333333334, "
        "43.999583333333334, 6.000416666666667, 45.000416666666666]"
    )
