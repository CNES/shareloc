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
This module contains functions to generate stereo-rectification epipolar grids
"""

# Standard imports
import math
from typing import Tuple, Union

# Third party imports
import numpy as np
import rasterio
from affine import Affine

# Shareloc imports
from shareloc.geofunctions.dtm_intersection import DTMIntersection
from shareloc.geofunctions.localization import Localization, coloc
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.image import Image
from shareloc.proj_utils import transform_index_to_physical_point


def write_epipolar_grid(grid, filename, xy_convention=True):
    """
    Write epipolar grid in a tiff file

    :param grid: epipolar grid
    :type grid: shareloc.Image
    :param filename:  output filename
    :type filename: string
    :param xy_convention: True: write grid with xy convention : [band 1 = col displacement, band 2 = row displacement]
        False: write grid with yx convention : [band 1 = row displacement, band 2 = col displacement]
    :param xy_convention: bool
    """
    band, row, col = grid.data.shape

    with rasterio.open(
        filename, "w", driver="GTiff", dtype=np.float64, width=col, height=row, count=band, transform=grid.transform
    ) as source_ds:
        if xy_convention:
            source_ds.write(grid.data[1, :, :], 1)
            source_ds.write(grid.data[0, :, :], 2)
        else:
            source_ds.write(grid.data[0, :, :], 1)
            source_ds.write(grid.data[1, :, :], 2)


def compute_epipolar_angle(end_line, start_line):
    """
    Define the epipolar angle

    :param end_line: ending of the epipolar line (georeferenced coordinates)
    :type end_line: 1D np.array [row, col, altitude] or 2D np.array [number of points, [row, col, altitude]]
    :param start_line: beginning of the epipolar line (georeferenced coordinates)
    :type start_line: 1D np.array [row, col, altitude] or 2D np.array [number of points, [row, col, altitude]]
    :return: epipolar angle
    :rtype: float or 1D np.array
    """
    # Only one point, expand the shape of the array
    if len(end_line.shape) == 1:
        end_line = np.expand_dims(end_line, axis=0)
        start_line = np.expand_dims(start_line, axis=0)

    # Compute the equation of the epipolar line y = a*x + b and define the epipolare angle
    alpha = np.zeros(end_line.shape[0])

    # Same columns, positive direction
    same_col_positive = (end_line[:, 1] == start_line[:, 1]) & (end_line[:, 0] > start_line[:, 0])
    alpha[same_col_positive] = 0.5 * math.pi

    # Same columns, negative direction
    same_col_negative = (end_line[:, 1] == start_line[:, 1]) & (end_line[:, 0] <= start_line[:, 0])
    alpha[same_col_negative] = -0.5 * math.pi

    # Different columns, positive direction
    diff_col_pos = np.where((end_line[:, 1] != start_line[:, 1]) & (end_line[:, 1] > start_line[:, 1]))
    slope = (end_line[diff_col_pos[0], 0] - start_line[diff_col_pos[0], 0]) / (
        end_line[diff_col_pos[0], 1] - start_line[diff_col_pos[0], 1]
    )
    alpha[diff_col_pos] = np.arctan(slope)

    # Different columns, negative direction
    diff_col_neg = np.where((end_line[:, 1] != start_line[:, 1]) & (end_line[:, 1] <= start_line[:, 1]))
    slope = (end_line[diff_col_neg[0], 0] - start_line[diff_col_neg[0], 0]) / (
        end_line[diff_col_neg[0], 1] - start_line[diff_col_neg[0], 1]
    )
    alpha[diff_col_neg[0]] = math.pi + np.arctan(slope)

    return alpha


def compute_local_epipolar_line(geom_model_left, geom_model_right, left_point, elevation, elevation_offset):
    """
    Estimate the beginning and the ending of local epipolar line in left image

    :param geom_model_left: geometric model of the left image
    :type geom_model_left: GeomodelTemplate
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: GeomodelTemplate
    :param left_point: georeferenced coordinates in the left image
    :type left_point: 1D numpy array : [row coord, col coord, altitude]
                      or 2D numpy array : (number of points, [row coord, col coord, altitude])
    :param elevation: elevation
    :type elevation: shareloc.dtm or float
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: float
    :return: Coordinates of the beginning and the ending of local epipolar line in the left image
    :rtype: Tuple(1D np.array [row, col, altitude], 1D numpy array [row, col, altitude])
            or Tuple(2D np.array (nb points, [row, col, altitude]), 2D np.array (nb points, [row, col, altitude]))
    """
    # Only one point, expand the shape of the array
    if len(left_point.shape) == 1:
        left_point = np.expand_dims(left_point, axis=0)

    # Right correspondent of the left coordinates
    right_corr = np.zeros((left_point.shape[0], 3))
    right_corr[:, 0], right_corr[:, 1], right_corr[:, 2] = coloc(
        geom_model_left, geom_model_right, left_point[:, 0], left_point[:, 1], elevation
    )
    ground_elev = np.array(right_corr[:, 2])

    # Find the beginning of the epipolar line in the left image, using right correspondent at lower elevation
    right_corr[:, 2] = ground_elev - elevation_offset
    epi_line_start = np.zeros((left_point.shape[0], 3))
    epi_line_start[:, 0], epi_line_start[:, 1], epi_line_start[:, 2] = coloc(
        geom_model_right, geom_model_left, right_corr[:, 0], right_corr[:, 1], right_corr[:, 2]
    )

    # Find the ending of the epipolar line in the left image, using right correspondent at higher elevation
    right_corr[:, 2] = ground_elev + elevation_offset
    epi_line_end = np.zeros((left_point.shape[0], 3))
    epi_line_end[:, 0], epi_line_end[:, 1], epi_line_end[:, 2] = coloc(
        geom_model_right, geom_model_left, right_corr[:, 0], right_corr[:, 1], right_corr[:, 2]
    )

    return epi_line_start, epi_line_end


# pylint: disable=too-many-locals
def prepare_rectification(left_im, geom_model_left, geom_model_right, elevation, epi_step, elevation_offset):
    """
    Determine size and spacing of the epipolar grids.
    Determine size of the epipolar images and the upper-left origin of the stereo-rectified left image (starting point)

    :param left_im: left image
    :type left_im: shareloc.image object
    :param geom_model_left: geometric model of the left image
    :type geom_model_left: GeomodelTemplate
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: GeomodelTemplate
    :param elevation: elevation
    :type elevation: shareloc.dtm or float
    :param epi_step: epipolar step
    :type epi_step: int
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: int
    :return:
        - epipolar grids spacing (pixel size), 1D np.array [row pixel size, col pixel size]
        - epipolar grids size, 1D np.array [number of row, number of columns]
        - epipolar images size, 1D np.array [number of row, number of columns]
        - epipolar grid corners in left image geometry [ul, ll, lr, ur]
            2D np.array [georef corner_row, georef corner_col, altitude]
    :rtype: Tuple
    """
    # Choose a square spacing
    mean_spacing = 0.5 * (abs(left_im.pixel_size_col) + abs(left_im.pixel_size_row))

    #  Pixel size (spacing) of the epipolar grids, convention [row, col]
    grid_pixel_size = np.full(shape=2, fill_value=epi_step, dtype=np.float64)
    grid_pixel_size *= mean_spacing

    # Georeferenced coordinates of the upper-left origin of left image : [row, col, altitude]
    origin_row, origin_col = transform_index_to_physical_point(left_im.transform, 0, 0)
    left_origin = np.array([origin_row, origin_col])

    # --- Compute the parameters of the local epipolar line at the left image origin ---
    local_epi_start, local_epi_end = compute_local_epipolar_line(
        geom_model_left, geom_model_right, left_origin, elevation, elevation_offset
    )
    local_epi_start = np.squeeze(local_epi_start)
    local_epi_end = np.squeeze(local_epi_end)

    # 2) Compute epipolar angle using the begin and the end of the left local epipolar line
    alpha = compute_epipolar_angle(local_epi_end, local_epi_start)[0]

    # 3) Compute unitary vectors, useful for moving in rows and columns in epipolar geometry
    # Unit vector tangent to the epipolar line (moving along line)
    unit_vector_along_epi_x = np.cos(alpha)
    unit_vector_along_epi_y = np.sin(alpha)
    # Unit vector orthogonal to epipolar direction (moving to next line)
    unit_vector_ortho_epi_x = -np.sin(alpha)
    unit_vector_ortho_epi_y = np.cos(alpha)

    # 4) Compute the bounding box of the left input image in the epipolar geometry
    # Coordinates of the 4 corners
    ulx = 0
    uly = 0
    urx = unit_vector_along_epi_x * left_im.nb_columns * left_im.pixel_size_col
    ury = unit_vector_ortho_epi_x * left_im.nb_columns * left_im.pixel_size_col
    llx = unit_vector_along_epi_y * left_im.nb_rows * left_im.pixel_size_row
    lly = unit_vector_ortho_epi_y * left_im.nb_rows * left_im.pixel_size_row
    lrx = (
        unit_vector_along_epi_x * left_im.nb_columns * left_im.pixel_size_col
        + unit_vector_along_epi_y * left_im.nb_rows * left_im.pixel_size_row
    )
    lry = (
        unit_vector_ortho_epi_x * left_im.nb_columns * left_im.pixel_size_col
        + unit_vector_ortho_epi_y * left_im.nb_rows * left_im.pixel_size_row
    )

    # Bounding box
    minx = min(urx, llx, lrx, ulx)
    miny = min(ury, lly, lry, uly)
    maxx = max(urx, llx, lrx, ulx)
    maxy = max(ury, lly, lry, uly)

    # 5) Compute the size of epipolar images
    rectified_image_size = [int((maxy - miny) / mean_spacing), int((maxx - minx) / mean_spacing)]

    # 6) Georeferenced coordinates of the [ul, ll, lr, ur] position of left epipolar image
    left_epi_ul = [
        left_origin[0] + (unit_vector_along_epi_y * minx + unit_vector_ortho_epi_y * miny),
        left_origin[1] + (unit_vector_along_epi_x * minx + unit_vector_ortho_epi_x * miny),
        (local_epi_start[2] + local_epi_end[2]) / 2.0,
    ]
    left_epi_lr = [
        left_origin[0] + (unit_vector_along_epi_y * (maxx + epi_step) + unit_vector_ortho_epi_y * (maxy + epi_step)),
        left_origin[1] + (unit_vector_along_epi_x * (maxx + epi_step) + unit_vector_ortho_epi_x * (maxy + epi_step)),
        (local_epi_start[2] + local_epi_end[2]) / 2.0,
    ]
    left_epi_ur = [
        left_origin[0] + (unit_vector_along_epi_y * minx + unit_vector_ortho_epi_y * (maxy + epi_step)),
        left_origin[1] + (unit_vector_along_epi_x * minx + unit_vector_ortho_epi_x * (maxy + epi_step)),
        (local_epi_start[2] + local_epi_end[2]) / 2.0,
    ]
    left_epi_ll = [
        left_origin[0] + (unit_vector_along_epi_y * (maxx + epi_step) + unit_vector_ortho_epi_y * miny),
        left_origin[1] + (unit_vector_along_epi_x * (maxx + epi_step) + unit_vector_ortho_epi_x * miny),
        (local_epi_start[2] + local_epi_end[2]) / 2.0,
    ]

    footprint = np.array([left_epi_ul, left_epi_ll, left_epi_lr, left_epi_ur])

    # 7) Compute the size of the epipolar grids, convention [nb_row, nb_col]
    # Two cells are added to the grid in order to harmonize the OTB conventions.
    grid_size = [int(rectified_image_size[0] / epi_step + 2), int(rectified_image_size[1] / epi_step + 2)]

    return grid_pixel_size, grid_size, rectified_image_size, footprint


def get_epipolar_extent(
    left_im, geom_model_left, geom_model_right, elevation=0.0, epi_step=30.0, elevation_offset=50.0, margin=0.0
):
    """
    return epipolar footprint using reprojection of epipolar geometry in left image.

    :param left_im: left image
    :type left_im: shareloc.image object
    :param geom_model_left: geometric model of the left image
    :type geom_model_left: GeomodelTemplate
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: GeomodelTemplate
    :param elevation: elevation
    :type elevation: shareloc.dtm or float
    :param epi_step: epipolar step
    :type epi_step: float
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: float
    :param margin: footprint margin (in degrees)
    :type margin: float
    :return: [lon_min,lat_min,lon max,lat max] (2D np.array)
    :rtype: numpy.array
    """
    __, __, __, footprint = prepare_rectification(
        left_im, geom_model_left, geom_model_right, elevation, epi_step, elevation_offset
    )

    loc_left = Localization(geom_model_left, image=left_im)
    footprint = footprint[:, 0:2]
    on_ground_pos = loc_left.direct(footprint[:, 0], footprint[:, 1], 0, using_geotransform=False)
    [lon_min, lat_min, __] = np.min(on_ground_pos, 0)
    [lon_max, lat_max, __] = np.max(on_ground_pos, 0)
    # Sometimes a margin is added because we don't know the epipolar grid footprint size.
    return np.array([lat_min - margin, lon_min - margin, lat_max + margin, lon_max + margin])


def initialize_grids(epi_step, nb_row, nb_col):
    """
    Initialize left and right epipolar grids : set geo-transform and zeros data

    :param epi_step: epipolar step
    :param nb_row: rows of the grid
    :param nb_col: columns of the grid
    :return: left epipolar grid, right epipolar grid
    :rtype: Tuple(shareloc.Image, shareloc.Image)
    """
    # Initialize left and right epipolar grids
    left_grid = Image(image_path=None)

    # Convention :
    # | col pixel size, row rotation , origin col upper-left|
    # | col rotation,   row pixel size,  , origin row upper-left|
    left_grid_geo_transform = np.array(
        [epi_step, 0, -(epi_step * 0.5), 0, epi_step, -(epi_step * 0.5)], dtype=np.float64
    )
    left_grid.set_metadata(nb_row, nb_col, 2, left_grid_geo_transform, datatype=np.float64)

    right_grid = Image(image_path=None)
    right_grid_geo_transform = np.array(
        [epi_step, 0, -(epi_step * 0.5), 0, epi_step, -(epi_step * 0.5)], dtype=np.float64
    )
    right_grid.set_metadata(nb_row, nb_col, 2, right_grid_geo_transform, datatype=np.float64)

    return left_grid, right_grid


def moving_along_axis(
    geom_model_left: GeoModelTemplate,
    geom_model_right: GeoModelTemplate,
    current_coords: np.ndarray,
    spacing: float,
    elevation: Union[float, DTMIntersection],
    epi_step: int,
    epi_angles: np.ndarray,
    axis: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Moving to the next line in epipolar geometry

    :param geom_model_left: geometric model of the left image
    :type geom_model_left: GeoModelTemplate
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: GeoModelTemplate
    :param current_coords: current line in the left epipolar geometry
        or current georeferenced coordinates in left epipolar line
    :type current_coords: 1D np.array [row, col, altitude]
        or 2D numpy array (number rows in epipolar geometry, [row, col, altitude])
    :param spacing: image spacing along axis dimension (in general mean image spacing)
    :type spacing: float
    :param elevation: elevation
    :type elevation: shareloc.dtm or float
    :param epi_step: epipolar step
    :type epi_step: int
    :param epi_angles: epipolar angle
    :type epi_angles: np.ndarray
    :param axis: displacement direction (0 = along columns, 1 = along lines)
    :type axis: int
    :return: left and right positions in epipolar grid
    :rtype: Tuple([row, col, altitude], [row, col, altitude])
        or Tuple(2D numpy array (number rows in epipolar geometry, [row, col, altitude]),
        2D numpy array (number rows in epipolar geometry, [row, col, altitude]))
    """

    if axis not in [0, 1]:
        raise ValueError(f"axis value {axis} is not available")

    epi_angles = epi_angles + (1 - axis) * np.pi / 2

    unit_vector_along_epi_x = epi_step * spacing * np.cos(epi_angles)
    unit_vector_along_epi_y = epi_step * spacing * np.sin(epi_angles)

    next_left = np.copy(current_coords)
    next_left[:, 0] += unit_vector_along_epi_y
    next_left[:, 1] += unit_vector_along_epi_x

    # Find the corresponding next pixels in the right image
    next_right = np.zeros(next_left.shape, dtype=next_left.dtype)
    next_right[:, 0], next_right[:, 1], next_right[:, 2] = coloc(
        geom_model_left, geom_model_right, next_left[:, 0], next_left[:, 1], elevation
    )

    return next_left, next_right


# disable for api symmetry between left and right data
# pylint: disable=unused-argument
# pylint: disable=too-many-locals
def compute_stereorectification_epipolar_grids(
    left_im, geom_model_left, right_im, geom_model_right, elevation=0.0, epi_step=1, elevation_offset=50.0
):
    """
    Compute stereo-rectification epipolar grids

    :param left_im: left image
    :type left_im: shareloc.image object
    :param geom_model_left: geometric model of the left image
    :type geom_model_left: GeomodelTemplate
    :param right_im: right image
    :type right_im: shareloc.image object
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: GeomodelTemplate
    :param elevation: elevation
    :type elevation: shareloc.dtm or float
    :param epi_step: epipolar step
    :type epi_step: int
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: float
    :return:
        - left epipolar grid, shareloc.image object convention [[row displacement, col displacement], nb rows, nb cols]
        - right epipolar grid, shareloc.image object convention [[row displacement, col displacement], nb rows, nb cols]
        - number of rows of the epipolar image, int
        - number of columns of the epipolar image, int
        - mean value of the baseline to sensor altitude ratio, float
    :rtype: Tuple
    """
    # Retrieve grids : spacing (pixel size) and size
    # Retrieve epipolar image : size and upper-left origin in the left image geometry (starting point)
    __, grid_size, rectified_image_size, footprint = prepare_rectification(
        left_im, geom_model_left, geom_model_right, elevation, epi_step, elevation_offset
    )

    # Use the mean spacing as before
    mean_spacing = 0.5 * (abs(left_im.pixel_size_col) + abs(left_im.pixel_size_row))

    left_grid, right_grid = initialize_grids(epi_step, grid_size[0], grid_size[1])

    # Starting points are the upper-left origin of the left epipolar image, and it's correspondent in the right image
    start_left = np.array(np.copy(footprint[0]))
    start_left = np.reshape(start_left, (1, -1))
    start_right = np.zeros(3, dtype=start_left.dtype)
    start_right = np.reshape(start_right, (1, -1))
    init_row, init_col, init_alt = coloc(
        geom_model_left, geom_model_right, start_left[:, 0], start_left[:, 1], elevation
    )
    # Convert ndarray coloc output into float 64 (Bug python3.9 et 3.10 not allowed anymore)
    # TODO: clean epipolar grids generation conversion globally with refacto/optimization
    start_right[:, 0] = init_row[0]
    start_right[:, 1] = init_col[0]
    start_right[:, 2] = init_alt[0]

    mean_baseline_ratio = 0

    # Compute the starting point of each epipolar line to be able to move along the lines (useful to vectorize the code)
    # Georeferenced coordinates of each starting epipolar lines in left and right image
    left_epi_lines = [np.copy(start_left)]
    right_epi_lines = [np.copy(start_right)]

    # For each rows of the epipolar geometry, define left and right starting coordinates of each epipolar lines
    for __ in range(grid_size[0] - 1):
        # --- Compute left local epipolar line, useful for moving to the next line ---
        local_epi_start, local_epi_end = compute_local_epipolar_line(
            geom_model_left, geom_model_right, left_epi_lines[-1], elevation, elevation_offset
        )

        # epipolar angle using the begin and the end of the left local epipolar line
        alpha = compute_epipolar_angle(local_epi_end, local_epi_start)
        # Find the start of next line in epipolar geometry
        next_epi_line_left, next_epi_line_right = moving_along_axis(
            geom_model_left,
            geom_model_right,
            left_epi_lines[-1],
            mean_spacing,
            elevation,
            epi_step,
            alpha,
            0,
        )

        # Save the starting points, useful to be able to move along the lines in the next loop
        left_epi_lines.append(np.copy(next_epi_line_left))
        right_epi_lines.append(np.copy(next_epi_line_right))

    # Left and right epipolar coordinates of the current point
    left_epi_coords = np.squeeze(np.array(left_epi_lines))
    right_epi_coords = np.squeeze(np.array(right_epi_lines))

    # Moving along epipolar lines
    rows = np.arange(grid_size[0])
    for col in range(grid_size[1]):
        # Estimate the displacement values of the current pixels
        # Cast grid index to georeferenced grid coordinates
        current_left_grid = transform_index_to_physical_point(left_grid.transform, rows, np.repeat(col, rows.shape[0]))
        current_right_grid = transform_index_to_physical_point(
            right_grid.transform, rows, np.repeat(col, rows.shape[0])
        )

        left_grid.data[0, :, col] = left_epi_coords[:, 0] - current_left_grid[0]
        left_grid.data[1, :, col] = left_epi_coords[:, 1] - current_left_grid[1]
        right_grid.data[0, :, col] = right_epi_coords[:, 0] - current_right_grid[0]
        right_grid.data[1, :, col] = right_epi_coords[:, 1] - current_right_grid[1]

        # Compute left local epipolar line, useful to estimate the local baseline ratio and moving to the next pixels
        local_epi_start, local_epi_end = compute_local_epipolar_line(
            geom_model_left, geom_model_right, left_epi_coords, elevation, elevation_offset
        )
        # Estimate the local baseline ratio
        local_baseline_ratio = np.sqrt(
            (local_epi_end[:, 1] - local_epi_start[:, 1]) * (local_epi_end[:, 1] - local_epi_start[:, 1])
            + (local_epi_end[:, 0] - local_epi_start[:, 0]) * (local_epi_end[:, 0] - local_epi_start[:, 0])
        ) / (2 * elevation_offset)
        mean_baseline_ratio += np.sum(local_baseline_ratio)

        # epipolar angle using the begin and the end of the left local epipolar lines
        alphas = compute_epipolar_angle(local_epi_end, local_epi_start)

        # Move to the next pixels in the epipolar line (moving along lines)
        left_epi_coords, right_epi_coords = moving_along_axis(
            geom_model_left, geom_model_right, left_epi_coords, mean_spacing, elevation, epi_step, alphas, 1
        )

    # Compute the mean baseline ratio
    mean_baseline_ratio /= grid_size[0] * grid_size[1]

    return left_grid, right_grid, rectified_image_size[0], rectified_image_size[1], mean_baseline_ratio


# pylint: disable=too-many-arguments
def compute_strip_of_epipolar_grid(
    geom_model_left: GeoModelTemplate,
    geom_model_right: GeoModelTemplate,
    left_positions_point: np.ndarray,
    right_positions_point: np.ndarray,
    spacing: float,
    axis: int,
    strip_size: int,
    epi_step: int = 1,
    elevation: Union[float, DTMIntersection] = 0.0,
    elevation_offset: float = 50.0,
    epipolar_angles: np.ndarray = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Compute stereo-rectification epipolar grids by strip. We start with a starting positions_point,
    and optional already computed epipolar_angles at these points, and we move strip_size times along axis direction
    (0 = along columns, 1 = along lines) to compute a strip of epipolar grid. positions_points is a (rows,cols,3)
    array containing (rows,cols) 3D coordinates (row,col,alt).  If axis==0 (resp. 1)
    rows (resp. cols) dimension must be 1.

    :param geom_model_left: geometric model of the left image
    :type geom_model_left: GeomodelTemplate
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: GeomodelTemplate
    :param left_positions_point: array of size (rows,cols,3) containing positions for left points
    :type left_positions_point: np.ndarray
    :param right_positions_point: array of size (rows,cols,3) containing positions for left points
    :type right_positions_point: np.ndarray
    :param spacing: image spacing along axis dimension (in general mean image spacing)
    :type spacing: float
    :param axis: displacement direction (0 = along columns, 1 = along lines)
    :type axis: int
    :param strip_size: desired size of grid along one axis for strip generation
    :type strip_size: int
    :param epi_step: epipolar grid sampling step
    :type epi_step: int
    :param elevation: elevation
    :type elevation: shareloc.dtm or float
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: float
    :param epipolar_angles: 2D array (rows,cols) containing epipolar angles at each grid node (angle in between epipolar
        coordinate system and left image coordinate system), size must be coherent with positions_point 1st,2nd dim
    :type epipolar_angles: np.ndarray
    :return:
        - left  epipolar positions grid in shape (rows,cols,3) rows (resp. cols) is strip_size if axis = 0 (resp. 1)
        - right epipolar positions grid in shape (rows,cols,3) rows (resp. cols) is strip_size if axis = 0 (resp. 1)
        - array of epipolar angles (rows,cols)
        - mean baseline ratio of the strip computed on rows x cols elements minus provided epipolar angles shape.
            We consider that if epipolar angles are provided, then baseline ratio for these elements have been
            already computed.
    :rtype: Tuple
    """

    # Instantiate mean baseline ratio
    mean_baseline_ratio = 0

    # compute the output size depending on axis
    size_shape = (
        (strip_size, left_positions_point.shape[1], 3) if axis == 0 else (left_positions_point.shape[0], strip_size, 3)
    )
    epi_angles_out = np.zeros((size_shape[0], size_shape[1]))

    # rename
    current_left_point = left_positions_point
    current_right_point = right_positions_point

    # copy first elements in output grids
    left_grid, right_grid = np.zeros(size_shape), np.zeros(size_shape)
    left_grid[0 : current_left_point.shape[0], 0 : current_left_point.shape[1], :] = current_left_point
    right_grid[0 : current_left_point.shape[0], 0 : current_left_point.shape[1], :] = current_right_point

    current_left_point = np.squeeze(current_left_point)

    # squeeze reduce 2 dimension if shape is like (1,1,3) we want to keep (1,3) dims
    if current_left_point.ndim == 1:
        current_left_point = current_left_point[np.newaxis, :]

    # if epipolar angles is not given as input we first compute the epipolar direction, otherwise
    # we consider that baseline ration has been already computed thus already_computed_ratio index is updated
    if epipolar_angles is None:
        already_computed_ratio = 0.0
        local_epi_start, local_epi_end = compute_local_epipolar_line(
            geom_model_left, geom_model_right, current_left_point, elevation, elevation_offset
        )

        epipolar_angles = compute_epipolar_angle(local_epi_end, local_epi_start)

        local_baseline_ratio = np.sqrt(
            (local_epi_end[:, 1] - local_epi_start[:, 1]) * (local_epi_end[:, 1] - local_epi_start[:, 1])
            + (local_epi_end[:, 0] - local_epi_start[:, 0]) * (local_epi_end[:, 0] - local_epi_start[:, 0])
        ) / (2 * elevation_offset)
        mean_baseline_ratio += np.sum(local_baseline_ratio)
    else:
        already_computed_ratio = epipolar_angles.shape[0] * epipolar_angles.shape[1]
    epipolar_angles = np.squeeze(epipolar_angles)

    if axis == 0:
        epi_angles_out[0, :] = epipolar_angles
    else:
        epi_angles_out[:, 0] = epipolar_angles

    # Grid creation
    # 1/ move along axis
    # 2/ compute local epipolar line
    # 3/ compute locale epipolar angle
    # 4/ fill the output
    # 5/ compute and increment the baseline ratio
    for point in range(1, strip_size):
        current_left_point, current_right_point = moving_along_axis(
            geom_model_left, geom_model_right, current_left_point, spacing, elevation, epi_step, epipolar_angles, axis
        )

        local_epi_start, local_epi_end = compute_local_epipolar_line(
            geom_model_left, geom_model_right, current_left_point, elevation, elevation_offset
        )

        epipolar_angles = compute_epipolar_angle(local_epi_end, local_epi_start)

        # Stock values in good shape
        if axis == 0:
            epi_angles_out[point, :] = epipolar_angles
            left_grid[point, :, :] = current_left_point
            right_grid[point, :, :] = current_right_point
        else:
            epi_angles_out[:, point] = epipolar_angles
            left_grid[:, point, :] = current_left_point
            right_grid[:, point, :] = current_right_point

        local_baseline_ratio = np.sqrt(
            (local_epi_end[:, 1] - local_epi_start[:, 1]) * (local_epi_end[:, 1] - local_epi_start[:, 1])
            + (local_epi_end[:, 0] - local_epi_start[:, 0]) * (local_epi_end[:, 0] - local_epi_start[:, 0])
        ) / (2 * elevation_offset)
        mean_baseline_ratio += np.sum(local_baseline_ratio)

    # Compute the mean baseline ratio
    mean_baseline_ratio /= left_grid.shape[0] * left_grid.shape[1] - already_computed_ratio
    return left_grid, right_grid, epi_angles_out, mean_baseline_ratio


def positions_to_displacement_grid(
    left_grid: np.ndarray, right_grid: np.ndarray, epi_step: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    :param left_grid: left epipolar positions grids
    :type left_grid: np.ndarray
    :param right_grid: right epipolar positions grids
    :type right_grid: np.ndarray
    :param epi_step: epipolar step
    :type epi_step: float
    :return:
        - left epipolar displacement grid, np.ndarray
        - right epipolar displacement grid, np.ndarray
    :rtype: Tuple(np.ndarray, np.ndarray)
    """
    rows = np.arange(left_grid.shape[0], dtype=float)
    cols = np.arange(left_grid.shape[1], dtype=float)

    cols_v, rows_v = np.meshgrid(cols, rows)
    transform = Affine(epi_step, 0, -(epi_step * 0.5), 0, epi_step, -(epi_step * 0.5))
    grids_cols, grids_rows = transform_index_to_physical_point(transform, cols_v, rows_v)

    left_grid[:, :, 0] = left_grid[:, :, 0] - grids_rows
    left_grid[:, :, 1] = left_grid[:, :, 1] - grids_cols
    right_grid[:, :, 0] = right_grid[:, :, 0] - grids_rows
    right_grid[:, :, 1] = right_grid[:, :, 1] - grids_cols

    return left_grid, right_grid
