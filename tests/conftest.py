#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).
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
Test module for rectification grid interpolation class shareloc/geofunctions/rectification*.py
Ground truth references (gt_{left/right}_grid*.tif) have been generated using OTB StereoRectificationGridGenerator
application.
"""
# pylint: disable=redefined-outer-name
# Standard imports
import os

# Third party imports
import pytest

# Shareloc imports
from shareloc.geofunctions.rectification import init_inputs_rectification
from shareloc.geomodels.geomodel import GeoModel
from shareloc.image import Image

# Shareloc test imports
from tests.helpers import data_path


@pytest.fixture()
def init_rpc_geom_model():
    """
    init geomodel fixture
    """
    geom_model_left = GeoModel(os.path.join(data_path(), "rectification", "left_image.geom"))
    geom_model_right = GeoModel(os.path.join(data_path(), "rectification", "right_image.geom"))

    return geom_model_left, geom_model_right


@pytest.fixture()
def init_inputs_rectification_fixture(init_rpc_geom_model):
    """
    Inputs rectification fixture
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0

    (
        left_position_point,
        right_position_point,
        spacing,
        grid_size,
        __,
    ) = init_inputs_rectification(
        left_im, geom_model_left, right_im, geom_model_right, default_elev, epi_step, elevation_offset
    )

    return left_position_point, right_position_point, spacing, epi_step, elevation_offset, default_elev, grid_size
