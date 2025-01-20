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
Test module for rectification grid interpolation class shareloc/geofunctions/rectification*.py
Ground truth references (otb_{left/right}_grid*.tif) have been generated using OTB StereoRectificationGridGenerator
application.
"""
# pylint: disable=duplicate-code
# Standard imports
import math
import os

# Third party imports
import numpy as np
import pytest
import rasterio

from shareloc.dtm_reader import dtm_reader

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection
from shareloc.geofunctions.rectification import (  # write_epipolar_grid,
    compute_epipolar_angle,
    compute_stereorectification_epipolar_grids,
    get_epipolar_extent,
    moving_along_axis,
    positions_to_displacement_grid,
    prepare_rectification,
)
from shareloc.geofunctions.rectification_grid import RectificationGrid
from shareloc.geomodels import GeoModel
from shareloc.image import Image

# Shareloc test imports
from tests.helpers import data_path


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_rpc(init_rpc_geom_model):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: default to 0.0
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, default_elev, epi_step, elevation_offset
    )
    # to deplacement grid
    left_grid, right_grid, _ = positions_to_displacement_grid(left_grid, right_grid, epi_step)

    # OTB reference
    otb_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid.tif")).read()
    otb_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid.tif")).read()

    # update baseline
    # write_epipolar_grid(left_grid, os.path.join(data_path(),'shareloc_gt_left_grid.tif'),grid_geotransform)
    # write_epipolar_grid(right_grid, os.path.join(data_path(),'shareloc_gt_right_grid.tif'), grid_geotransform)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert otb_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1e-2)
    assert otb_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1e-2)

    assert otb_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1e-2)
    assert otb_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1e-2)

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.704004705
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)

    # Shareloc non regression reference
    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "shareloc_gt_left_grid.tif")).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_right_grid.tif")
    ).read()

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1], right_grid[:, :, 0])
    np.testing.assert_allclose(reference_right_grid[0], right_grid[:, :, 1], atol=2.0e-12)

    # Check mean_baseline_ratio
    reference_mean_br = 0.7040047235162911
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-11)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_rpc_displacement(init_rpc_geom_model):
    """
    Test epipolar grids generation as displacement grids :
     check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: default to 0.0
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im,
        geom_model_left,
        right_im,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        as_displacement_grid=True,
    )

    # OTB reference
    otb_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid.tif")).read()
    otb_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid.tif")).read()

    # update baseline
    # write_epipolar_grid(left_grid, os.path.join(data_path(),'shareloc_gt_left_grid.tif'),grid_geotransform)
    # write_epipolar_grid(right_grid, os.path.join(data_path(),'shareloc_gt_right_grid.tif'), grid_geotransform)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert otb_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1e-2)
    assert otb_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1e-2)

    assert otb_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1e-2)
    assert otb_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1e-2)

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.704004705
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)

    # Shareloc non regression reference
    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "shareloc_gt_left_grid.tif")).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_right_grid.tif")
    ).read()

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1], right_grid[:, :, 0])
    np.testing.assert_allclose(reference_right_grid[0], right_grid[:, :, 1], atol=2.0e-12)

    # Check mean_baseline_ratio
    reference_mean_br = 0.7040047235162911
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-11)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_rpc_with_margins(init_rpc_geom_model):
    """
    Test epipolar grids generation a with margin :
     check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: default to 0.0
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    margin = 5
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im,
        geom_model_left,
        right_im,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        margin=margin,
    )

    # update baseline
    # origin = [float(-margin), float(-margin)]
    # spacing = [float(epi_step), float(epi_step)]
    # grid_geotransform = (
    #     origin[0] - 0.5 * spacing[0],
    #     spacing[0],
    #     0.0,
    #     origin[1] - 0.5 * spacing[1],
    #     0.0,
    #     spacing[1],
    # )
    # write_epipolar_grid(
    #     left_grid,
    #     os.path.join(data_path(), "rectification", 'shareloc_gt_left_grid_with_margins.tif'),
    #     grid_geotransform
    # )
    # write_epipolar_grid(
    #     right_grid,
    #     os.path.join(data_path(), "rectification", 'shareloc_gt_right_grid_with_margins.tif'),
    #     grid_geotransform
    # )

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.704004673
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)

    # Shareloc non regression reference
    reference_left_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_left_grid_with_margins.tif")
    ).read()
    reference_right_grid = rasterio.open(
        os.path.join(data_path(), "rectification", "shareloc_gt_right_grid_with_margins.tif")
    ).read()

    # Check epipolar grids
    np.testing.assert_array_equal(reference_left_grid[1], left_grid[:, :, 0])
    np.testing.assert_array_equal(reference_left_grid[0], left_grid[:, :, 1])

    np.testing.assert_array_equal(reference_right_grid[1], right_grid[:, :, 0])
    np.testing.assert_allclose(reference_right_grid[0], right_grid[:, :, 1], atol=2.0e-12)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_rpc_alti(init_rpc_geom_model):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: alti=100
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 100.0
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, default_elev, epi_step, elevation_offset
    )
    # to deplacement grid
    left_grid, right_grid, _ = positions_to_displacement_grid(left_grid, right_grid, epi_step)

    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid_100.tif")).read()
    reference_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid_100.tif")).read()

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert reference_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1e-2)
    assert reference_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1e-2)

    assert reference_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1e-2)
    assert reference_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1e-2)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.7039927244
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_rpc_dtm_geoid(init_rpc_geom_model):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation: SRTM DTM + Geoid egm96_15
    """

    # first instantiate geometric models left and right (here RPC geometrics model)
    geom_model_left, geom_model_right = init_rpc_geom_model

    # read the images
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    # we use DTM and Geoid: a DTMIntersection class has to be used
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(dtm_file, geoid_file)
    dtm_ventoux = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    # compute rectification grid sampled at 30 pixels
    epi_step = 30
    elevation_offset = 50
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, dtm_ventoux, epi_step, elevation_offset
    )
    # to deplacement grid
    left_grid, right_grid, _ = positions_to_displacement_grid(left_grid, right_grid, epi_step)

    # evaluate the results by comparison with OTB
    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid_dtm.tif")).read()
    reference_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid_dtm.tif")).read()

    # baseline update if necessary
    # write_epipolar_grid(left_grid, os.path.join(data_path(),'grid_left_dtm.tif'))
    # write_epipolar_grid(right_grid, os.path.join(data_path(),'grid_right_dtm.tif'))

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert reference_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1e-2)
    assert reference_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1e-2)

    assert reference_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1e-2)
    assert reference_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1e-2)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.7039416432
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_rpc_dtm_geoid_roi(init_rpc_geom_model):
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: RPC
    Earth elevation ROI: SRTM DTM + Geoid egm96_15
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    extent = get_epipolar_extent(left_im, geom_model_left, geom_model_right, additional_margin=0.0016667)
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=extent,
        roi_is_in_physical_space=True,
        fill_nodata=None,
        fill_value=0.0,
    )
    dtm_ventoux = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    epi_step = 30
    elevation_offset = 50
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, dtm_ventoux, epi_step, elevation_offset
    )

    # to deplacement grid
    left_grid, right_grid, _ = positions_to_displacement_grid(left_grid, right_grid, epi_step)

    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid_dtm.tif")).read()
    reference_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid_dtm.tif")).read()

    # update baseline
    # write_epipolar_grid(left_grid, os.path.join(data_path(),'grid_left_dtm_roi.tif'))
    # write_epipolar_grid(right_grid, os.path.join(data_path(),'grid_right_dtm_roi.tif'))

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]

    assert reference_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1e-2)
    assert reference_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1e-2)
    assert reference_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1e-2)
    assert reference_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1e-2)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.7039416432
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-5)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_grid():
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio


    Input Geomodels: Grid
    Earth elevation: default to 0.0
    """

    # first instantiate geometric models left and right (here Grid geometric model)
    geom_model_left = GeoModel(
        os.path.join(data_path(), "grid/phr_ventoux/GRID_PHR1B_P_201308051042194_SEN_690908101-001.tif"), "GRID"
    )
    geom_model_right = GeoModel(
        os.path.join(data_path(), "grid/phr_ventoux/GRID_PHR1B_P_201308051042523_SEN_690908101-002.tif"), "GRID"
    )

    # read the images
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    default_elev = 0.0

    # compute rectification grid sampled at 30 pixels
    epi_step = 30
    elevation_offset = 50
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, default_elev, epi_step, elevation_offset
    )

    # to deplacement grid
    left_grid, right_grid, _ = positions_to_displacement_grid(left_grid, right_grid, epi_step)

    # OTB reference
    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid.tif")).read()
    reference_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid.tif")).read()

    # update baseline
    # write_epipolar_grid(left_grid, os.path.join(data_path(),'grid_left_elev_0.tif'))
    # write_epipolar_grid(right_grid, os.path.join(data_path(),'grid_right_elev_0.tif'))

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert reference_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1.2e-2)
    assert reference_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1.2e-2)

    assert reference_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1.2e-2)
    assert reference_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1.2e-2)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    referecne_mean_br = 0.704004705
    assert mean_br == pytest.approx(referecne_mean_br, abs=1e-4)


@pytest.mark.unit_tests
def test_compute_stereorectification_epipolar_grids_geomodel_grid_dtm_geoid():
    """
    Test epipolar grids generation : check epipolar grids, epipolar image size, mean_baseline_ratio

    Input Geomodels: Grids
    Earth elevation: SRTM DTM + Geoid egm96_15
    """

    # first instantiate geometric models left and right (here Grid geometrics model)
    geom_model_left = GeoModel(
        os.path.join(data_path(), "grid/phr_ventoux/GRID_PHR1B_P_201308051042194_SEN_690908101-001.tif"), "GRID"
    )
    geom_model_right = GeoModel(
        os.path.join(data_path(), "grid/phr_ventoux/GRID_PHR1B_P_201308051042523_SEN_690908101-002.tif"), "GRID"
    )

    # read the images
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    # we use DTM and Geoid a DTMIntersection class has to be used
    dtm_file = os.path.join(data_path(), "dtm", "srtm_ventoux", "srtm90_non_void_filled", "N44E005.hgt")
    geoid_file = os.path.join(data_path(), "dtm", "geoid", "egm96_15.gtx")
    dtm_image = dtm_reader(
        dtm_file,
        geoid_file,
        roi=None,
        roi_is_in_physical_space=True,
        fill_nodata="mean",
        fill_value=0.0,
    )
    dtm_ventoux = DTMIntersection(
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )

    # compute rectification grid sampled at 30 pixels
    epi_step = 30
    elevation_offset = 50
    left_grid, right_grid, [img_size_row, img_size_col], mean_br = compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, dtm_ventoux, epi_step, elevation_offset
    )

    # to deplacement grid
    left_grid, right_grid, _ = positions_to_displacement_grid(left_grid, right_grid, epi_step)

    # evaluate the results by comparison with OTB
    reference_left_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid_dtm.tif")).read()
    reference_right_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_right_grid_dtm.tif")).read()

    # baseline update if necessary
    # write_epipolar_grid(left_grid, os.path.join(data_path(),'grid_left_dtm.tif'))
    # write_epipolar_grid(right_grid, os.path.join(data_path(),'grid_right_dtm.tif'))

    # Check epipolar grids
    # OTB convention is [col, row], shareloc convention is [row, col]
    assert reference_left_grid[1] == pytest.approx(left_grid[:, :, 0], abs=1e-2)
    assert reference_left_grid[0] == pytest.approx(left_grid[:, :, 1], abs=1e-2)

    assert reference_right_grid[1] == pytest.approx(right_grid[:, :, 0], abs=1e-2)
    assert reference_right_grid[0] == pytest.approx(right_grid[:, :, 1], abs=1e-2)

    # Check size of rectified images
    assert img_size_row == 612
    assert img_size_col == 612

    # Check mean_baseline_ratio
    # ground truth mean baseline ratio from OTB
    reference_mean_br = 0.7039416432
    assert mean_br == pytest.approx(reference_mean_br, abs=1e-4)


@pytest.mark.parametrize("row,col", [(15, 0)])
@pytest.mark.unit_tests
def test_rectification_grid_interpolation_one_point(row, col):
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(data_path(), "rectification_grids", f"grid_{id_scene_right}.tif")
    rectif_grid = RectificationGrid(grid_filename, is_displacement_grid=True)
    # value at position [15,15]
    value_row = np.sum(rectif_grid.row_positions[0, 0:2]) / 2.0
    value_col = np.sum(rectif_grid.col_positions[0, 0:2]) / 2.0
    coords = rectif_grid.interpolate(np.array((col, row)))
    assert value_row == pytest.approx(coords[0, 1], abs=1e-4)
    assert value_col == pytest.approx(coords[0, 0], abs=1e-4)


@pytest.mark.unit_tests
def test_rectification_grid_interpolation():
    """
    Test interpolation on rectification grid
    """
    id_scene_right = "P1BP--2017092838319324CP"
    grid_filename = os.path.join(data_path(), "rectification_grids", f"grid_{id_scene_right}.tif")

    rectif_grid = RectificationGrid(grid_filename, is_displacement_grid=True)
    # value at position [15,15]

    value_row = np.sum(rectif_grid.row_positions[0:2, 0:2]) / 4.0
    value_col = np.sum(rectif_grid.col_positions[0:2, 0:2]) / 4.0
    sensor_positions = np.zeros((2, 2))
    sensor_positions[0, :] = [15.0, 15.0]
    sensor_positions[1, :] = [0.0, 0.0]
    coords = rectif_grid.interpolate(sensor_positions)
    assert value_col == pytest.approx(coords[0, 0], abs=1e-4)
    assert value_row == pytest.approx(coords[0, 1], abs=1e-4)


@pytest.mark.unit_tests
def test_rectification_grid_extrapolation():
    """
    Test interpolation on rectification grid
    """
    grid_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid_ventoux.tif")
    rectif_grid = RectificationGrid(grid_filename, is_displacement_grid=True)

    sensor_positions = np.zeros((2, 2))
    sensor_positions[0, :] = [30.0, -10.0]
    sensor_positions[1, :] = [30.0, 691.0]
    coords = rectif_grid.interpolate(sensor_positions)
    assert pytest.approx(coords[0, 1], abs=1e-10) == 5065.72347005208303016
    assert pytest.approx(coords[1, 1], abs=1e-10) == 4883.84894205729142413474619389


@pytest.mark.unit_tests
def test_prepare_rectification(init_rpc_geom_model):
    """
    Test prepare rectification : check grids size, epipolar image size, and left epipolar starting point
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    margin = 0
    __, grid_size, rectified_image_size, footprint, _ = prepare_rectification(
        left_im,
        geom_model_left,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        margin,
    )

    # check size of the epipolar grids
    assert grid_size[0] == 22
    assert grid_size[1] == 22

    # check size of rectified images
    assert rectified_image_size[0] == 612
    assert rectified_image_size[1] == 612

    # check the first epipolar point in the left image
    # ground truth values from OTB
    otb_output_origin_in_left_image = [5625.78139593690139008685946465, 5034.15635707952696975553408265, 0]

    # OTB convention is [col, row, altitude], shareloc convention is [row, col, altitude]
    assert otb_output_origin_in_left_image[1] == pytest.approx(footprint[0][0], abs=1e-5)
    assert otb_output_origin_in_left_image[0] == pytest.approx(footprint[0][1], abs=1e-5)

    assert footprint[0][2] == otb_output_origin_in_left_image[2]


@pytest.mark.unit_tests
def test_prepare_rectification_footprint(init_rpc_geom_model):
    """
    Test prepare rectification : check footprint
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    margin = 0
    _, _, _, footprint, _ = prepare_rectification(
        left_im,
        geom_model_left,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        margin,
    )

    ground_truth = np.array(
        [
            [5034.15635485, 5625.78139208, 0.0],
            [5654.75411307, 5459.06023724, 0.0],
            [5488.03295823, 4838.46247902, 0.0],
            [4867.43520001, 5005.18363386, 0.0],
        ]
    )

    assert np.all(ground_truth == pytest.approx(footprint, abs=1e-5))


@pytest.mark.unit_tests
def test_prepare_rectification_footprint_with_margins(init_rpc_geom_model):
    """
    Test prepare rectification : check footprint
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    margin = 5
    _, _, _, footprint, _ = prepare_rectification(
        left_im,
        geom_model_left,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        margin,
    )

    ground_truth = np.array(
        [
            [4928.20978948, 5809.56203656, 0.0],
            [5838.53475755, 5565.00680260, 0.0],
            [5593.97952359, 4654.68183453, 0.0],
            [4683.65455552, 4899.23706849, 0.0],
        ]
    )

    assert np.all(ground_truth == pytest.approx(footprint, abs=1e-5))


@pytest.mark.unit_tests
def test_rectification_moving_along_lines(init_rpc_geom_model):
    """
    Test moving along line in epipolar geometry
    """

    geom_model_left, geom_model_right = init_rpc_geom_model

    current_coords = np.array([[5000.5, 5000.5, 0.0]], dtype=np.float64)
    mean_spacing = 1
    epi_step = 1
    alphas = 0
    default_elev = 0.0
    # ground truth next pixel
    # col pixel size of the image
    col_pixel_size = 1.0
    reference_next_cords = np.array([[5000.5, 5000.5 + col_pixel_size, 0.0]], dtype=np.float64)

    next_cords, _ = moving_along_axis(
        geom_model_left,
        geom_model_right,
        current_coords,
        mean_spacing,
        default_elev,
        epi_step,
        alphas,
        1,
    )

    np.testing.assert_array_equal(reference_next_cords, next_cords)


@pytest.mark.unit_tests
def test_rectification_moving_to_next_line(init_rpc_geom_model):
    """
    Test moving to next line in epipolar geometry
    """

    geom_model_left, geom_model_right = init_rpc_geom_model

    current_coords = np.array([[5000.5, 5000.5, 0.0]], dtype=np.float64)
    mean_spacing = 1
    epi_step = 1
    alphas = 0
    default_elev = 0.0
    # ground truth next pixel
    # row pixel size of the image
    row_pixel_size = 1.0
    reference_next_cords = np.array([[5000.5 + row_pixel_size, 5000.5, 0.0]], dtype=np.float64)

    next_cords, _ = moving_along_axis(
        geom_model_left,
        geom_model_right,
        current_coords,
        mean_spacing,
        default_elev,
        epi_step,
        alphas,
        0,
    )

    np.testing.assert_array_equal(reference_next_cords, next_cords)


@pytest.mark.unit_tests
def test_rectification_moving_to_axis_error():
    """
    Test moving along axis with wrong axis
    """
    geom_model_left = GeoModel(os.path.join(data_path(), "rectification", "left_image.geom"))
    geom_model_right = GeoModel(os.path.join(data_path(), "rectification", "right_image.geom"))

    current_coords = np.array([5000.5, 5000.5, 0.0], dtype=np.float64)
    mean_spacing = 1
    epi_step = 1
    alphas = 0
    default_elev = 0.0

    with pytest.raises(ValueError):
        _, _ = moving_along_axis(
            geom_model_left, geom_model_right, current_coords, mean_spacing, default_elev, epi_step, alphas, 2
        )


@pytest.mark.unit_tests
def test_epipolar_angle():
    """
    test epipolar angle computation
    """
    # First case : same column, positive direction [row, col, alt]
    start_line_1 = np.array([1, 0, 0])
    end_line_1 = np.array([2, 0, 0])

    reference_alpha_1 = math.pi / 2.0
    alpha = compute_epipolar_angle(end_line_1, start_line_1)
    assert alpha == reference_alpha_1

    # Second case : same column, negative direction [row, col, alt]
    start_line_2 = np.array([2, 0, 0])
    end_line_2 = np.array([1, 0, 0])

    reference_alpha_2 = -(math.pi / 2.0)
    alpha = compute_epipolar_angle(end_line_2, start_line_2)
    assert alpha == reference_alpha_2

    # Third case : different column, positive direction [row, col, alt]
    start_line_3 = np.array([2, 0, 0])
    end_line_3 = np.array([1, 1, 0])

    slope = (1 - 2) / (1 - 0)

    reference_alpha_3 = np.arctan(slope)
    alpha = compute_epipolar_angle(end_line_3, start_line_3)
    assert alpha == reference_alpha_3

    # Fourth case : different column, negative direction [row, col, alt]
    start_line_4 = np.array([2, 1, 0])
    end_line_4 = np.array([1, 0, 0])

    slope = (1 - 2) / (0 - 1)
    reference_alpha_4 = math.pi + np.arctan(slope)
    alpha = compute_epipolar_angle(end_line_4, start_line_4)
    assert alpha == reference_alpha_4

    # With multiple point
    start_lines = np.stack((start_line_1, start_line_2, start_line_3, start_line_4))
    end_lines = np.stack((end_line_1, end_line_2, end_line_3, end_line_4))
    reference_alphas = np.stack((reference_alpha_1, reference_alpha_2, reference_alpha_3, reference_alpha_4))

    alphas = compute_epipolar_angle(end_lines, start_lines)
    np.testing.assert_array_equal(alphas, reference_alphas)


@pytest.mark.unit_tests
def test_rectification_grid_pos_inside_prepare_footprint_bounding_box(init_rpc_geom_model):
    """
    Test that epipolar grid is inside the footprint returned by prepare_rectification
    """
    # Generate epipolar grid and parameters
    # Ground truth generated by GridBasedResampling function from OTB.
    epi_grid = rasterio.open(os.path.join(data_path(), "rectification", "otb_left_grid.tif"))
    transform = epi_grid.transform
    # Grid origin
    origin_row = transform[5]
    origin_col = transform[2]
    # Pixel size
    pixel_size_row = transform[4]
    pixel_size_col = transform[0]
    # Get array
    grid_col_dep = epi_grid.read(1)
    grid_row_dep = epi_grid.read(2)
    # Origin coordinates
    grid_origin_row = origin_row + pixel_size_row / 2.0
    grid_origin_col = origin_col + pixel_size_col / 2.0
    # Create grid with displacements
    grid_pos_col = np.arange(grid_origin_col, grid_origin_col + epi_grid.width * pixel_size_col, step=pixel_size_col)
    grid_pos_col = np.tile(grid_pos_col, (epi_grid.height, 1))
    grid_pos_row = np.arange(grid_origin_row, grid_origin_row + epi_grid.height * pixel_size_row, step=pixel_size_row)
    grid_pos_row = np.tile(grid_pos_row, (epi_grid.width, 1)).T
    pos_col = grid_col_dep + grid_pos_col
    pos_row = grid_row_dep + grid_pos_row

    # Get positions
    positions = np.stack((pos_col.flatten(), pos_row.flatten()))  # X, Y

    # Get grid footprint
    grid_footprint = np.array(
        [
            positions[:, int(grid_origin_row)],
            positions[:, int(grid_origin_row + grid_pos_row.shape[0] - 1)],
            positions[:, int(grid_origin_row + grid_pos_row.size - 1)],
            positions[:, int(grid_origin_row + grid_pos_row.size - 22)],
        ]
    )
    # OTB convention is [col, row, altitude], shareloc convention is [row, col, altitude]
    grid_footprint = grid_footprint[:, [1, 0]]

    # Compute shareloc epipolar footprint
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0
    margin = 0
    _, _, _, footprint, _ = prepare_rectification(
        left_im,
        geom_model_left,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        margin,
    )
    footprint = footprint[:, 0:2]

    # Bounding box
    min_row, max_row = min(footprint[:, 0]), max(footprint[:, 0])
    min_col, max_col = min(footprint[:, 1]), max(footprint[:, 1])

    # Test that grid_footprint is in epipolar footprint
    assert np.all(np.logical_and(min_row < grid_footprint[:, 0], grid_footprint[:, 0] < max_row))
    assert np.all(np.logical_and(min_col < grid_footprint[:, 1], grid_footprint[:, 1] < max_col))
