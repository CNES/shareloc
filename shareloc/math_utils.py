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
This module contains the mathematical functions for shareloc
"""

import numpy as np


def interpol_bilin(mats, nb_rows, nb_cols, delta_shift_row, delta_shift_col):
    """
    bilinear interpolation on multi layer  matrix
    :param mats: multi layer grid (: , nb_rows,nb_cols)
    :param nb_rows: line number of mats
    :param nb_cols: column number of mats
    :param delta_shift_row: position (line)
    :param delta_shift_col: position (column)
    :return interpolated value on each layer
    :rtype list
    """
    if delta_shift_row < 0:
        lower_shift_row = 0
    elif delta_shift_row >= nb_rows - 1:
        lower_shift_row = nb_rows - 2
    else:
        lower_shift_row = int(np.floor(delta_shift_row))
    upper_shift_row = lower_shift_row + 1

    if delta_shift_col < 0:
        lower_shift_col = 0
    elif delta_shift_col >= nb_cols - 1:
        lower_shift_col = nb_cols - 2
    else:
        lower_shift_col = int(np.floor(delta_shift_col))
    upper_shift_col = lower_shift_col + 1
    # (col_shift, row_shift) are subpixel distance to interpolate along each axis
    col_shift = delta_shift_col - lower_shift_col
    row_shift = delta_shift_row - lower_shift_row
    # Altitude
    matis = []
    for mat in mats:
        mati = (
            (1 - col_shift) * (1 - row_shift) * mat[:, lower_shift_row, lower_shift_col]
            + col_shift * (1 - row_shift) * mat[:, lower_shift_row, upper_shift_col]
            + (1 - col_shift) * row_shift * mat[:, upper_shift_row, lower_shift_col]
            + col_shift * row_shift * mat[:, upper_shift_row, upper_shift_col]
        )
        matis.append(mati)
    return matis


def interpol_bilin_vectorized(mats, nb_rows, nb_cols, delta_shift_row, delta_shift_col):
    """
    bilinear interpolation on multi points and layer  matrix
    :param mats: multi layer grid (: , nb_rows,nb_cols)
    :type mats : list
    :param nb_rows: line number of mats
    :type nb_rows : int
    :param nb_cols: column number of mats
    :type nb_cols : int
    :param delta_shift_row: position (line)
    :type delta_shift_row: 1D numpy.ndarray, dtype=float64
    :param delta_shift_col: position (column)
    :type delta_shift_col: 1D numpy.ndarray, dtype=float64
    :return interpolated value on each layer
    :rtype list
    """
    lower_shift_row = np.floor(delta_shift_row).astype(int)
    lower_shift_row[delta_shift_row < 0] = 0
    lower_shift_row[delta_shift_row >= (nb_rows - 1)] = nb_rows - 2
    upper_shift_row = np.copy(lower_shift_row) + 1

    lower_shift_col = np.floor(delta_shift_col).astype(int)
    lower_shift_col[delta_shift_col < 0] = 0
    lower_shift_col[delta_shift_col >= (nb_cols - 1)] = nb_cols - 2
    upper_shift_col = np.copy(lower_shift_col) + 1

    # (col_shift,row_shift) are subpixel distance to interpolate along each axis
    col_shift = delta_shift_col - lower_shift_col
    row_shift = delta_shift_row - lower_shift_row

    # Altitude
    matis = []
    for mat in mats:
        mati = (
            (1 - col_shift) * (1 - row_shift) * mat[:, lower_shift_row, lower_shift_col]
            + col_shift * (1 - row_shift) * mat[:, lower_shift_row, upper_shift_col]
            + (1 - col_shift) * row_shift * mat[:, upper_shift_row, lower_shift_col]
            + col_shift * row_shift * mat[:, upper_shift_row, upper_shift_col]
        )
        matis.append(mati)

    return matis
