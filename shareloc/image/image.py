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
Image class to handle Image data.
"""

import rasterio
import numpy as np
from affine import Affine


class Image:
    """ class Image to handle image data """

    def __init__(self, image_path, read_data=False):
        """
        constructor
        :param image_path : image path
        :type image_path  : string or None
        :param read_data  : read image data
        :type read_data  : bool
        """
        if image_path is not None:
            # Image path
            self.image_path = image_path

            # Rasterio dataset
            self.dataset = rasterio.open(image_path)

            # Georeferenced coordinates of the upper-left origin
            self.origin_row = self.dataset.transform[5]
            self.origin_col = self.dataset.transform[2]

            # Image size
            self.nb_rows = self.dataset.height
            self.nb_columns = self.dataset.width

            # Pixel size
            self.pixel_size_row = self.dataset.transform[4]
            self.pixel_size_col = self.dataset.transform[0]

            # Geo-transform of type Affine with convention :
            # | pixel size col,   row rotation, origin col |
            # | col rotation  , pixel size row, origin row |
            self.transform = self.dataset.transform

            self.data = None
            if read_data:
                # Data of shape (nb band, nb row, nb col)
                self.data = self.dataset.read()

    def set_metadata(self, nb_row, nb_col, nb_band, transform, datatype=np.float32):
        """
        Set metadata and create data (np.array filled with zeros of shape (nb_band, nb_row, nb_col))

        :param nb_row: number of row
        :type nb_row: int
        :param nb_col: number of col
        :type nb_col: int
        :param nb_band: number of band
        :type nb_band: int
        :param transform: geo-transform
        :type transform: 1D numpy array with convention :
                | pixel size col, row rotation, origin col, col rotation, pixel size row, origin row |
        :param datatype: data type
        :type datatype: numpy type
        """
        self.transform = Affine(transform[0], transform[1], transform[2], transform[3], transform[4], transform[5])

        # Georeferenced coordinates of the upper-left origin
        self.origin_row = transform[5]
        self.origin_col = transform[2]

        # Image size
        self.nb_rows = nb_row
        self.nb_columns = nb_col

        # Pixel size
        self.pixel_size_row = transform[4]
        self.pixel_size_col = transform[0]

        self.data = np.zeros((nb_band, nb_row, nb_col), dtype=datatype)

    def transform_index_to_physical_point(self, row, col):
        """
        Transform index to physical point

        :param row: row index
        :type row: int or 1D numpy array
        :param col: col index
        :type col: int or 1D numpy array
        :return: Georeferenced coordinates (row, col)
        :rtype: Tuple(georeference row float or 1D np.array, georeference col float or 1D np.array)
        """

        # trans with convention : | pixel size col, row rotation, origin col, col rotation, pixel size row, origin row |
        trans = self.transform
        col_geo = trans[2] + (col + 0.5) * trans[0] + (row + 0.5) * trans[1]
        row_geo = trans[5] + (col + 0.5) * trans[3] + (row + 0.5) * trans[4]

        return row_geo, col_geo

    def transform_physical_point_to_index(self, row_geo, col_geo):
        """
        Transform physical point to index

        :param row_geo: physical point row
        :type row_geo: int or 1D numpy array
        :param col_geo: physical point col
        :type col_geo: int or 1D numpy array
        :return: index coordinates (row, col)
        :rtype: Tuple(row float or 1D np.array, col float or 1D np.array)
        """
        row = (row_geo - self.origin_row) / self.pixel_size_row - 0.5
        col = (col_geo - self.origin_col) / self.pixel_size_col - 0.5
        return row, col