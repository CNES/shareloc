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
Tests dedicated to rectification footprint coverage.
These tests verify that the rectification footprint covers both the left
and the right images, for constant elevation and DTM-based rectification.
"""
# pylint: disable=duplicate-code
# Standard imports
import os

# Third party imports
import numpy as np
import pytest

from shareloc.dtm_reader import dtm_reader

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection
from shareloc.geofunctions.localization import coloc
from shareloc.geofunctions.rectification import (  # write_epipolar_grid,
    _get_image_corners,
    _transform_points_to_epipolar_frame,
    compute_epipolar_angle,
    compute_local_epipolar_line,
    compute_stereorectification_epipolar_grids,
    get_epipolar_extent,
    prepare_rectification,
)
from shareloc.image import Image
from shareloc.proj_utils import transform_index_to_physical_point

# Shareloc test imports
from tests.helpers import data_path


def _compute_test_epipolar_frame(
    left_im,
    geom_model_left,
    geom_model_right,
    elevation,
    elevation_offset,
    footprint,
):
    """
    Reconstruct the local epipolar frame used by prepare_rectification.

    This helper converts points expressed in the left-image geometry into the
    rectified image coordinate system, allowing the coverage of the rectified
    extent to be verified.
    """
    mean_spacing = 0.5 * (abs(left_im.pixel_size_col) + abs(left_im.pixel_size_row))

    origin_row, origin_col = transform_index_to_physical_point(
        left_im.transform,
        0,
        0,
    )
    left_origin = np.array(
        [origin_row, origin_col],
        dtype=float,
    )

    local_epi_start, local_epi_end, _ = compute_local_epipolar_line(
        geom_model_left,
        geom_model_right,
        left_origin,
        elevation,
        elevation_offset,
    )

    local_epi_start = np.squeeze(local_epi_start)
    local_epi_end = np.squeeze(local_epi_end)

    alpha = compute_epipolar_angle(
        local_epi_end,
        local_epi_start,
    )[0]

    unit_vector_along_epi_x = np.cos(alpha)
    unit_vector_along_epi_y = np.sin(alpha)
    unit_vector_ortho_epi_x = -np.sin(alpha)
    unit_vector_ortho_epi_y = np.cos(alpha)

    footprint_epipolar = _transform_points_to_epipolar_frame(
        footprint,
        left_origin,
        unit_vector_along_epi_x,
        unit_vector_along_epi_y,
        unit_vector_ortho_epi_x,
        unit_vector_ortho_epi_y,
    )

    return {
        "origin": left_origin,
        "along_x": unit_vector_along_epi_x,
        "along_y": unit_vector_along_epi_y,
        "ortho_x": unit_vector_ortho_epi_x,
        "ortho_y": unit_vector_ortho_epi_y,
        "minx": float(np.min(footprint_epipolar[:, 0])),
        "miny": float(np.min(footprint_epipolar[:, 1])),
        "mean_spacing": mean_spacing,
    }


def _transform_points_to_rectified_frame(
    points,
    epipolar_frame,
):
    """
    Transform points from left-image physical coordinates to rectified
    image coordinates.

    :param points: points in left-image geometry, shape (N, 2) or (N, 3),
        convention [row, column, ...]
    :type points: np.ndarray
    :param epipolar_frame: local epipolar frame parameters returned by
        _compute_test_epipolar_frame
    :type epipolar_frame: dict
    :return: rectified coordinates, shape (N, 2), convention [row, column]
    :rtype: np.ndarray
    """
    points_epi = _transform_points_to_epipolar_frame(
        points,
        epipolar_frame["origin"],
        epipolar_frame["along_x"],
        epipolar_frame["along_y"],
        epipolar_frame["ortho_x"],
        epipolar_frame["ortho_y"],
    )

    rectified_columns = (points_epi[:, 0] - epipolar_frame["minx"]) / epipolar_frame["mean_spacing"]

    rectified_rows = (points_epi[:, 1] - epipolar_frame["miny"]) / epipolar_frame["mean_spacing"]

    return np.column_stack(
        (
            rectified_rows,
            rectified_columns,
        )
    )


def _is_point_in_convex_polygon(point, polygon, tolerance=1e-8):
    """
    Check whether a point lies inside a convex polygon.

    :param point: point coordinates [row, col]
    :type point: np.ndarray
    :param polygon: ordered polygon vertices, shape (N, 2)
    :type polygon: np.ndarray
    :param tolerance: numerical tolerance
    :type tolerance: float
    :return: True if the point is inside or on the polygon boundary
    :rtype: bool
    """
    cross_products = []

    for index, first_vertex in enumerate(polygon):
        second_vertex = polygon[(index + 1) % len(polygon)]

        edge = second_vertex - first_vertex
        point_vector = point - first_vertex

        cross_product = edge[0] * point_vector[1] - edge[1] * point_vector[0]
        cross_products.append(cross_product)

    cross_products = np.asarray(cross_products)

    return bool(np.all(cross_products >= -tolerance) or np.all(cross_products <= tolerance))


@pytest.mark.unit_tests
def test_prepare_rectification_covers_right_image(init_rpc_geom_model):
    """
    Check that the rectification footprint covers the right image
    projected into the left image geometry at constant elevation.
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    epi_step = 30
    elevation_offset = 50
    default_elev = 0.0

    _, _, _, footprint, _ = prepare_rectification(
        left_im,
        geom_model_left,
        right_im,
        geom_model_right,
        default_elev,
        epi_step,
        elevation_offset,
        margin=0,
    )

    right_corners = _get_image_corners(right_im)

    right_corners_in_left, _ = coloc(
        geom_model_right,
        geom_model_left,
        right_corners[:, 0],
        right_corners[:, 1],
        elevation=default_elev,
        image1=right_im,
        image2=left_im,
    )

    assert all(
        _is_point_in_convex_polygon(
            right_corner[:2],
            footprint[:, :2],
        )
        for right_corner in right_corners_in_left
    )


@pytest.mark.unit_tests
def test_prepare_rectification_covers_right_image_dtm(
    init_rpc_geom_model,
):
    """
    Check that the rectification footprint covers the right image
    projected into the left image geometry for the minimum and
    maximum DTM elevations.
    """
    left_im = Image(os.path.join(data_path(), "rectification", "left_image.tif"))
    right_im = Image(os.path.join(data_path(), "rectification", "right_image.tif"))

    geom_model_left, geom_model_right = init_rpc_geom_model

    dtm_file = os.path.join(
        data_path(),
        "dtm",
        "srtm_ventoux",
        "srtm90_non_void_filled",
        "N44E005.hgt",
    )
    geoid_file = os.path.join(
        data_path(),
        "dtm",
        "geoid",
        "egm96_15.gtx",
    )

    extent = get_epipolar_extent(
        left_im,
        geom_model_left,
        right_im,
        geom_model_right,
        additional_margin=0.0016667,
    )

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

    (
        _,
        _,
        rectified_image_size,
        _,
    ) = compute_stereorectification_epipolar_grids(
        left_im,
        geom_model_left,
        right_im,
        geom_model_right,
        dtm_ventoux,
        epi_step,
        elevation_offset,
    )

    img_size_row, img_size_col = rectified_image_size

    valid_altitudes = np.asarray(dtm_image.alt_data, dtype=float)
    valid_altitudes = valid_altitudes[np.isfinite(valid_altitudes)]

    assert valid_altitudes.size > 0

    minimum_altitude = float(np.min(valid_altitudes))
    maximum_altitude = float(np.max(valid_altitudes))

    (
        _,
        _,
        _,
        footprint,
        _,
    ) = prepare_rectification(
        left_im,
        geom_model_left,
        right_im,
        geom_model_right,
        dtm_ventoux,
        epi_step,
        elevation_offset,
        margin=0,
    )

    epipolar_frame = _compute_test_epipolar_frame(
        left_im,
        geom_model_left,
        geom_model_right,
        dtm_ventoux,
        elevation_offset,
        footprint,
    )

    right_corners = _get_image_corners(right_im)

    right_rows = right_corners[:, 0]
    right_columns = right_corners[:, 1]

    right_corners_at_min_alt, _ = coloc(
        geom_model_right,
        geom_model_left,
        right_rows,
        right_columns,
        elevation=np.full(
            right_rows.shape,
            minimum_altitude,
            dtype=float,
        ),
        image1=right_im,
        image2=left_im,
    )

    right_corners_at_max_alt, _ = coloc(
        geom_model_right,
        geom_model_left,
        right_rows,
        right_columns,
        elevation=np.full(
            right_rows.shape,
            maximum_altitude,
            dtype=float,
        ),
        image1=right_im,
        image2=left_im,
    )

    right_corners_min_rectified = _transform_points_to_rectified_frame(
        right_corners_at_min_alt,
        epipolar_frame,
    )

    right_corners_max_rectified = _transform_points_to_rectified_frame(
        right_corners_at_max_alt,
        epipolar_frame,
    )

    coverage_tolerance = 1e-6

    def is_inside_rectified_image(point):
        """Check whether a point lies inside the rectified-image extent."""
        row, column = point

        return (
            -coverage_tolerance <= row <= img_size_row + coverage_tolerance
            and -coverage_tolerance <= column <= img_size_col + coverage_tolerance
        )

    assert all(is_inside_rectified_image(point) for point in right_corners_min_rectified)

    assert all(is_inside_rectified_image(point) for point in right_corners_max_rectified)
