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
Image class to handle Image data.
Shareloc Reference raster image input based on rasterio.
"""

# Standard imports
import logging

# Third party imports
import numpy as np
import rasterio
from affine import Affine


# pylint: disable=too-many-instance-attributes
class Image:
    """
    class Image to handle image data
    Shareloc reference for raster image input based on rasterio.
    """

    def __init__(self, image_path, read_data=False, roi=None, roi_is_in_physical_space=False):
        """
        constructor

        :param image_path : image path
        :type image_path  : string or None
        :param read_data  : read image data
        :type read_data  : bool
        :param roi  : region of interest [row_min,col_min,row_max,col_max] or [ymin,xmin,ymax,xmax] if
             roi_is_in_physical_space activated
        :type roi  : list
        :param roi_is_in_physical_space  : roi value in physical space
        :type roi_is_in_physical_space  : bool
        """
        if image_path is not None:
            # Image path
            self.image_path = image_path

            # Rasterio dataset
            self.dataset = rasterio.open(image_path)

            # Geo-transform of type Affine with convention :
            # | pixel size col,   row rotation, origin col |
            # | col rotation  , pixel size row, origin row |
            self.transform = self.dataset.transform

            roi_window = None
            if roi is not None:
                if roi_is_in_physical_space:
                    row_off, col_off = self.transform_physical_point_to_index(roi[0], roi[1])
                    row_max, col_max = self.transform_physical_point_to_index(roi[2], roi[3])
                    # index is relative to pixel center, here the roi is defined with corners
                    row_off += 0.5
                    col_off += 0.5
                    row_max += 0.5
                    col_max += 0.5
                    # in case of negative pixel size y
                    if row_off > row_max:
                        row_max, row_off = row_off, row_max
                    row_off = np.floor(row_off)
                    col_off = np.floor(col_off)
                    width = int(np.ceil(col_max - col_off))
                    height = int(np.ceil(row_max - row_off))
                    logging.info("roi in image , offset : %s %s size %s %s", col_off, row_off, width, height)
                else:
                    row_off = roi[0]
                    col_off = roi[1]
                    width = roi[3] - roi[1]
                    height = roi[2] - roi[0]
                roi_window = rasterio.windows.Window(col_off, row_off, width, height)
                self.transform = self.dataset.window_transform(roi_window)
                self.nb_rows = height
                self.nb_columns = width
            else:
                # Image size
                self.nb_rows = self.dataset.height
                self.nb_columns = self.dataset.width
                roi_window = None

            # Georeferenced coordinates of the upper-left origin
            self.origin_row = self.transform[5]
            self.origin_col = self.transform[2]

            # Pixel size
            self.pixel_size_row = self.transform[4]
            self.pixel_size_col = self.transform[0]

            # row/col rotation
            self.pixel_rotation_row = self.transform[3]
            self.pixel_rotation_col = self.transform[1]

            if self.dataset.crs is not None:
                self.epsg = self.dataset.crs.to_epsg()
            else:
                self.epsg = None

            self.nodata = self.dataset.nodata

            self.mask = None
            self.data = None
            if read_data:
                # Data of shape (nb band, nb row, nb col)
                self.data = np.squeeze(self.dataset.read(window=roi_window))
                if self.nodata is not None:
                    self.mask = np.squeeze(self.dataset.read_masks(window=roi_window))
                    logging.info("image contains %d nodata values ", np.sum(self.mask == 0))

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
        :type row: float or 1D np.array
        :param col: col index
        :type col: float or 1D np.array
        :return: Georeferenced coordinates (row, col)
        :rtype: Tuple(georeference row float or 1D np.array, georeference col float or 1D np.array)
        """

        trans = self.transform
        col_geo, row_geo = trans * (col + 0.5, row + 0.5)
        return row_geo, col_geo

    def transform_physical_point_to_index(self, row_geo, col_geo):
        """
        Transform physical point to index

        :param row_geo: physical point row
        :type row_geo: float or 1D np.array
        :param col_geo: physical point col
        :type col_geo: float or 1D np.array
        :return: index coordinates (row, col)
        :rtype: Tuple(row float or 1D np.array, col float or 1D np.array)
        """
        trans_inv = ~self.transform
        col, row = trans_inv * (col_geo, row_geo)
        return row - 0.5, col - 0.5
