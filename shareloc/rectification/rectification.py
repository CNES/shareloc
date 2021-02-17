#!/usr/bin/env python

# coding: utf8
#
# Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
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

import math
import numpy as np
import rasterio

from shareloc.localization import coloc
from shareloc.image.image import Image


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
    :rtype : float or 1D np.array
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

    return np.squeeze(alpha)


def get_local_altitude(dem=None, default_elev=0.0):
    """
    Get local altitude of point

    :param dem: use dem
    :type dem: shareloc.dtm
    :param default_elev: default altitude above ellipsoid
    :type default_elev: float
    :return: the local altitude of the point
    :rtype: float
    """
    # Set the default altitude above ellipsoid in case no information is available
    local_altitude = default_elev

    # Compute the altitude above ellipsoid of a geographic point
    if dem is not None:
        print("dem handling not implemented")
        # pass

    return local_altitude


def compute_local_epipolar_line(geom_model_left, geom_model_right, left_point, dem, default_elev, elevation_offset):
    """
    Estimate the beginning and the ending of local epipolar line in left image

    :param geom_model_left: geometric model of the left image
    :type geom_model_left: shareloc.grid or shareloc.rpc
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: shareloc.grid or shareloc.rpc
    :param left_point: georeferenced coordinates in the left image
    :type left_point: 1D numpy array : [row coord, col coord, altitude]
                      or 2D numpy array : (number of points, [row coord, col coord, altitude])
    :param dem: use dem for localisation
    :type dem: shareloc.dtm
    :param default_elev: default elevation (if dem is not None default elev is ignored)
    :type default_elev: float
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: int
    :return: Coordinates of the beginning and the ending of local epipolar line in the left image
    :rtype: Tuple(1D np.array [row, col, altitude], 1D numpy array [row, col, altitude])
            or Tuple(2D np.array (nb points, [row, col, altitude]), 2D np.array (nb points, [row, col, altitude]))
    """
    # Only one point, expand the shape of the array
    if len(left_point.shape) == 1:
        left_point = np.expand_dims(left_point, axis=0)
    left_point[:, 2] = get_local_altitude(dem, default_elev)
    # Right correspondent of the left coordinates
    right_corr = np.zeros((left_point.shape[0], 3))
    right_corr[:, 0], right_corr[:, 1], right_corr[:, 2] = coloc(
        geom_model_left, geom_model_right, left_point[:, 0], left_point[:, 1], left_point[:, 2]
    )

    # Find the beginning of the epipolar line in the left image, using right correspondent at lower elevation
    right_corr[:, 2] = left_point[:, 2] - elevation_offset
    epi_line_start = np.zeros((left_point.shape[0], 3))
    epi_line_start[:, 0], epi_line_start[:, 1], epi_line_start[:, 2] = coloc(
        geom_model_right, geom_model_left, right_corr[:, 0], right_corr[:, 1], right_corr[:, 2]
    )

    # Find the ending of the epipolar line in the left image, using right correspondent at higher elevation
    right_corr[:, 2] = left_point[:, 2] + elevation_offset
    epi_line_end = np.zeros((left_point.shape[0], 3))
    epi_line_end[:, 0], epi_line_end[:, 1], epi_line_end[:, 2] = coloc(
        geom_model_right, geom_model_left, right_corr[:, 0], right_corr[:, 1], right_corr[:, 2]
    )

    return np.squeeze(epi_line_start), np.squeeze(epi_line_end)


def prepare_rectification(left_im, geom_model_left, geom_model_right, dem, default_elev, epi_step, elevation_offset):
    """
    Determine size and spacing of the epipolar grids.
    Determine size of the epipolare images and the upper-left origin of the stereo-rectified left image (starting point)

    :param left_im: left image
    :type left_im: shareloc.image object
    :param geom_model_left: geometric model of the left image
    :type geom_model_left: shareloc.grid or shareloc.rpc
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: shareloc.grid or shareloc.rpc
    :param dem: use dem for localisation
    :type dem: shareloc.dtm
    :param default_elev: default elevation (if dem is not None default elev is ignored)
    :type default_elev: float
    :param epi_step: epipolar step
    :type epi_step: int
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: int
    :return: return :
        - epipolar grids spacing (pixel size), 1D np.array [row pixel size, col pixel size]
        - epipolar grids size, 1D np.array [number of row, number of columns]
        - epipolare images size, 1D np.array [number of row, number of columns]
        - epipolare image starting point in the left image, 1D np.array [georef row, georef col, altitude]

    :rtype: Tuple
    """
    # Choose a square spacing
    mean_spacing = 0.5 * (abs(left_im.pixel_size_col) + abs(left_im.pixel_size_row))

    #  Pixel size (spacing) of the epipolar grids, convention [row, col]
    grid_pixel_size = np.full(shape=2, fill_value=epi_step, dtype=np.float64)
    grid_pixel_size *= mean_spacing

    local_altitude = get_local_altitude(dem, default_elev)

    # Georeferenced coordinates of the upper-left origin of left image : [row, col, altitude]
    origin_row, origin_col = left_im.transform_index_to_physical_point(0, 0)
    left_origin = np.array([origin_row, origin_col, local_altitude])

    # --- Compute the parameters of the local epipolar line at the left image origin ---
    local_epi_start, local_epi_end = compute_local_epipolar_line(
        geom_model_left, geom_model_right, left_origin, dem, default_elev, elevation_offset
    )

    # 2) Compute epipolar angle using the begin and the end of the left local epipolar line
    alpha = compute_epipolar_angle(local_epi_end, local_epi_start)

    # 3) Compute unitary vectors, useful for moving in rows and columns in epipolar geometry
    # Unit vector tangent to the epipolar line (moving along line)
    unit_vector_along_epi_x = math.cos(alpha)
    unit_vector_along_epi_y = math.sin(alpha)
    # Unit vector orthogonal to epipolar direction (moving to next line)
    unit_vector_ortho_epi_x = -math.sin(alpha)
    unit_vector_ortho_epi_y = math.cos(alpha)

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
    minx = min(min(min(urx, llx), lrx), ulx)
    miny = min(min(min(ury, lly), lry), uly)
    maxx = max(max(max(urx, llx), lrx), ulx)
    maxy = max(max(max(ury, lly), lry), uly)

    # 5) Compute the size of epipolar images
    rectified_image_size = [int((maxy - miny) / mean_spacing), int((maxx - minx) / mean_spacing)]

    # 6) Georeferenced coordinates of the upper-left origin of left epiolar image (starting point)
    left_epi_origin = [
        left_origin[0] + (unit_vector_along_epi_y * minx + unit_vector_ortho_epi_y * miny),
        left_origin[1] + (unit_vector_along_epi_x * minx + unit_vector_ortho_epi_x * miny),
        local_altitude,
    ]

    # 7) Compute the size of the epipolar grids, convention [nb_row, nb_col]
    grid_size = [int(rectified_image_size[0] / epi_step + 2), int(rectified_image_size[1] / epi_step + 2)]

    return grid_pixel_size, grid_size, rectified_image_size, left_epi_origin


def initialize_grids(epi_step, nb_row, nb_col):
    """
    Initialize left and right epipolar grids : set geo-transform and zeros data

    :param epi_step: epipolar step
    :param nb_row: rows of the grid
    :param nb_col: columns of the grid
    :return: left epipolar grid, right epipolar grid
    :rtype : Tuple(shareloc.Image, shareloc.Image)
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


def moving_to_next_line(
    geom_model_left, geom_model_right, current_line, mean_spacing, dem, default_elev, epi_step, alpha
):
    """
    Moving to the next line in epipolar geometry

    :param geom_model_left: geometric model of the left image
    :type geom_model_left: shareloc.grid or  shareloc.rpc
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: shareloc.grid or  shareloc.rpc
    :param current_line: current line in the left epipolar geometry
    :type current_line: 1D np.array [row, col, altitude]
    :param mean_spacing: mean spacing of epipolar grids
    :type mean_spacing: int
    :param dem: use dem for localisation
    :type dem: shareloc.dtm
    :param default_elev: default elevation (if dem is not None default elev is ignored)
    :type default_elev: float
    :param epi_step: epipolar step
    :type epi_step: int
    :param alpha: epipolar angle
    :type alpha: float
    :return: left and right coordinates of the next line in epioplar geometry
    :rtype: Tuple([row, col, altitude], [row, col, altitude])
    """
    # Find the start of next line in epipolar geometry
    # Unit vector orthogonal to epipolar direction
    unit_vector_ortho_epi_y = epi_step * mean_spacing * math.cos(alpha)
    unit_vector_ortho_epi_x = -1 * mean_spacing * epi_step * math.sin(alpha)

    next_line_start_left = np.zeros(3, dtype=np.float64)

    next_line_start_left[0] = np.copy(current_line[0]) + unit_vector_ortho_epi_y
    next_line_start_left[1] = np.copy(current_line[1]) + unit_vector_ortho_epi_x
    local_altitude = get_local_altitude(dem, default_elev)
    next_line_start_left[2] = local_altitude

    # Find the corresponding starting point in the right image
    next_line_start_right = np.zeros(3, dtype=np.float64)
    next_line_start_right[0], next_line_start_right[1], next_line_start_right[2] = coloc(
        geom_model_left, geom_model_right, next_line_start_left[0], next_line_start_left[1], next_line_start_left[2]
    )
    return next_line_start_left, next_line_start_right


def moving_along_lines(
    geom_model_left, geom_model_right, current_left_coords, mean_spacing, dem, default_elev, epi_step, alphas
):
    """
    Determine the position of next pixels in the epipolare lines

    :param geom_model_left: geometric model of the left image
    :type geom_model_left: shareloc.grid or  shareloc.rpc
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: shareloc.grid or  shareloc.rpc
    :param current_left_coords: current georeferenced coordinates in left epipolare line
    :type current_left_coords: 2D numpy array (number rows in epipolar geometry, [row, col, altitude])
    :param mean_spacing: mean spacing of epipolar grids
    :type mean_spacing: int
    :param dem: use dem for localisation
    :type dem: shareloc.dtm
    :param default_elev: default elevation (if dem is not None default elev is ignored)
    :type default_elev: float
    :param epi_step: epipolar step
    :type epi_step: int
    :param alphas: epipolar angles of each local epipolar lines
    :type alphas: List
    :return: coordinates of next pixels in left epipolar line, coordinates of next pixels in right epipolar line
    :rtype: Tuple(2D numpy array (number rows in epipolar geometry, [row, col, altitude]),
                 2D numpy array (number rows in epipolar geometry, [row, col, altitude]))
    """
    # Determine position of next pixels in the epipolar line of the left image (moving along lines)
    # Unit vector tangent to the epipolar line
    unit_vector_along_epi_y = epi_step * mean_spacing * np.sin(alphas)
    unit_vector_along_epi_x = epi_step * mean_spacing * np.cos(alphas)

    # Move to the next pixels in left image
    next_left_coords = np.copy(current_left_coords)
    next_left_coords[:, 0] += unit_vector_along_epi_y
    next_left_coords[:, 1] += unit_vector_along_epi_x
    local_altitude = get_local_altitude(dem, default_elev)
    next_left_coords[:, 2] = local_altitude

    # Find the corresponding next pixels in the right image
    next_right_coords = np.zeros(next_left_coords.shape, dtype=next_left_coords.dtype)
    next_right_coords[:, 0], next_right_coords[:, 1], next_right_coords[:, 2] = coloc(
        geom_model_left, geom_model_right, next_left_coords[:, 0], next_left_coords[:, 1], next_left_coords[:, 2]
    )
    return next_left_coords, next_right_coords


def compute_stereorectification_epipolar_grids(
    left_im, geom_model_left, right_im, geom_model_right, dem=None, default_elev=0.0, epi_step=1, elevation_offset=50.0
):
    """
    Compute stereo-rectification epipolar grids

    :param left_im: left image
    :type left_im: shareloc.image object
    :param geom_model_left: geometric model of the left image
    :type geom_model_left: shareloc.grid or  shareloc.rpc
    :param right_im: right image
    :type right_im: shareloc.image object
    :param geom_model_right: geometric model of the right image
    :type geom_model_right: shareloc.grid or  shareloc.rpc
    :param dem: use dem for localisation
    :type dem: shareloc.dtm
    :param default_elev: default elevation (if dem is not None default elev is ignored)
    :type default_elev: float
    :param epi_step: epipolar step
    :type epi_step: int
    :param elevation_offset: elevation difference used to estimate the local tangent
    :type elevation_offset: float
    :return: return :
        - left epipolar grid, shareloc.image object convention [[row displacement, col displacement], nb rows, nb cols]
        - right epipolar grid, shareloc.image object convention [[row displacement, col displacement], nb rows, nb cols]
        - number of rows of the epipolar image, int
        - number of columns of the epipolar image, int
        - mean value of the baseline to sensor altitude ratio, float
    :rtype: Tuple
    """
    # Retrieve grids : spacing (pixel size) and size
    # Retrieve epipolar image : size and upper-left origin in the left image geometry (starting point)
    __, grid_size, rectified_image_size, left_epi_origin = prepare_rectification(
        left_im, geom_model_left, geom_model_right, dem, default_elev, epi_step, elevation_offset
    )

    # Use the mean spacing as before
    mean_spacing = 0.5 * (abs(left_im.pixel_size_col) + abs(left_im.pixel_size_row))

    left_grid, right_grid = initialize_grids(epi_step, grid_size[0], grid_size[1])

    # Starting points are the upper-left origin of the left epipolar image, and it's correspondent in the right image
    start_left = np.copy(left_epi_origin)
    start_left[2] = get_local_altitude(dem, default_elev)
    start_right = np.zeros(3, dtype=start_left.dtype)
    start_right[0], start_right[1], start_right[2] = coloc(
        geom_model_left, geom_model_right, start_left[0], start_left[1], start_left[2]
    )

    mean_baseline_ratio = 0

    # Compute the starting point of each epipolar line to be able to move along the lines (useful to vectorize the code)
    # Georeferenced coordinates of each starting epipolar lines in left and right image
    left_epi_lines = [np.copy(start_left)]
    right_epi_lines = [np.copy(start_right)]

    # For each rows of the epipolar geometry, define left and right starting coordinates of each epipolar lines
    for __ in range(grid_size[0] - 1):
        # --- Compute left local epipolar line, useful for moving to the next line ---
        local_epi_start, local_epi_end = compute_local_epipolar_line(
            geom_model_left, geom_model_right, left_epi_lines[-1], dem, default_elev, elevation_offset
        )

        # epipolar angle using the begin and the end of the left local epipolar line
        alpha = compute_epipolar_angle(local_epi_end, local_epi_start)
        # Find the start of next line in epipolar geometry
        next_epi_line_left, next_epi_line_right = moving_to_next_line(
            geom_model_left, geom_model_right, left_epi_lines[-1], mean_spacing, dem, default_elev, epi_step, alpha
        )

        # Save the starting points, useful to be able to move along the lines in the next loop
        left_epi_lines.append(np.copy(next_epi_line_left))
        right_epi_lines.append(np.copy(next_epi_line_right))

    # Left and right epipolar coordinates of the current point
    left_epi_coords = np.array(left_epi_lines)
    right_epi_coords = np.array(right_epi_lines)

    # Moving along epipolar lines
    rows = np.arange(grid_size[0])
    for col in range(grid_size[1]):
        # Estimate the displacement values of the current pixels
        # Cast grid index to georeferenced grid coordinates
        current_left_grid = left_grid.transform_index_to_physical_point(rows, np.repeat(col, rows.shape[0]))
        current_right_grid = right_grid.transform_index_to_physical_point(rows, np.repeat(col, rows.shape[0]))

        left_grid.data[0, :, col] = left_epi_coords[:, 0] - current_left_grid[0]
        left_grid.data[1, :, col] = left_epi_coords[:, 1] - current_left_grid[1]
        right_grid.data[0, :, col] = right_epi_coords[:, 0] - current_right_grid[0]
        right_grid.data[1, :, col] = right_epi_coords[:, 1] - current_right_grid[1]

        # Compute left local epipolar line, useful to estimate the local baseline ratio and moving to the next pixels
        local_epi_start, local_epi_end = compute_local_epipolar_line(
            geom_model_left, geom_model_right, left_epi_coords, dem, default_elev, elevation_offset
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
        left_epi_coords, right_epi_coords = moving_along_lines(
            geom_model_left, geom_model_right, left_epi_coords, mean_spacing, dem, default_elev, epi_step, alphas
        )

    # Compute the mean baseline ratio
    mean_baseline_ratio /= grid_size[0] * grid_size[1]

    return left_grid, right_grid, rectified_image_size[0], rectified_image_size[1], mean_baseline_ratio
