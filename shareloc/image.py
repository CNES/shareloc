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
# pylint: disable=no-member

# Standard imports
import logging

# Third party imports
import numpy as np
import rasterio
from affine import Affine

# Project import
from shareloc.proj_utils import transform_physical_point_to_index


# pylint: disable=too-many-instance-attributes
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-branches
class Image:
    """
    class Image to handle image data
    Shareloc reference for raster image input based on rasterio.
    """

    def __init__(
        self, image_path, read_data=False, roi=None, roi_is_in_physical_space=False, vertical_direction="auto"
    ):
        """
        constructor

        :param image_path : image path
        :type image_path  : string or None
        :param read_data  : read image data
        :type read_data  : bool
        :param roi  : region of interest [row_min,col_min,row_max,col_max] or [ymin,xmin,ymax,xmax] if
             roi_is_in_physical_space activated
        :type roi  : list
        :param roi_is_in_physical_space  : ROI value in physical space
        :type roi_is_in_physical_space  : bool
        :vertical_direction: option to choose direction when moving to next row, by default ("auto") there is
            no forced direction
                "north": downward, y>0
                "south": upward, y<0
        :type vertical_direction: str
        """

        self.vertical_direction = vertical_direction

        if image_path is not None:
            # Image path
            self.image_path = image_path

            # Rasterio dataset
            self.dataset = rasterio.open(image_path)

            # Geo-transform of type Affine with convention :
            # | pixel size col,   row rotation, origin col |
            # | col rotation  , pixel size row, origin row |
            self.transform = self.dataset.transform
            # bitwise not inversion (Affine.__invert implemented, pylint bug)
            self.trans_inv = ~self.transform  # pylint: disable=invalid-unary-operand-type
            if roi is not None:
                # User have set ROI in physical space or not
                if roi_is_in_physical_space:
                    row_off, col_off = transform_physical_point_to_index(self.trans_inv, roi[0], roi[1])
                    row_max, col_max = transform_physical_point_to_index(self.trans_inv, roi[2], roi[3])
                    # index is relative to pixel center, here the roi is defined with corners
                    row_off += 0.5
                    col_off += 0.5
                    row_max += 0.5
                    col_max += 0.5
                    # in case of negative pixel size y
                    if row_off > row_max:
                        row_max, row_off = row_off, row_max

                    row_off = max(np.floor(row_off), 0)
                    col_off = max(np.floor(col_off), 0)
                    width = int(np.ceil(col_max - col_off))
                    height = int(np.ceil(row_max - row_off))

                    logging.debug("roi in image , offset : %s %s size %s %s", col_off, row_off, width, height)
                else:
                    row_off = max(roi[0], 0)
                    col_off = max(roi[1], 0)
                    row_off = min(row_off, self.dataset.height)
                    col_off = min(col_off, self.dataset.width)
                    width = roi[3] - col_off
                    height = roi[2] - row_off

                # set boundaries for ROI
                width = int(min(width, self.dataset.width - col_off))
                height = int(min(height, self.dataset.height - row_off))

                roi_window = rasterio.windows.Window(col_off, row_off, width, height)
                self.transform = self.dataset.window_transform(roi_window)
                # bitwise not inversion (Affine.__invert implemented, pylint bug)
                self.trans_inv = ~self.transform  # pylint: disable=invalid-unary-operand-type
                self.nb_rows = height
                self.nb_columns = width
            else:
                # Image size if no ROI
                self.nb_rows = self.dataset.height
                self.nb_columns = self.dataset.width
                roi_window = None

            # Invert y axis if needed
            if self.vertical_direction == "north" and self.transform[4] < 0:  # force y positive
                logging.info("Changing y (vertical) axis direction: north y>0")
                self.transform = Affine(
                    self.transform[0],
                    self.transform[1],
                    self.transform[2],
                    -self.transform[3],
                    -self.transform[4],
                    -self.transform[5],
                )
                self.trans_inv = ~self.transform  # pylint: disable=invalid-unary-operand-type

            elif self.vertical_direction == "south" and self.transform[4] > 0:  # force y negative
                logging.info("Changing y (vertical) axis direction: south y<0")
                self.transform = Affine(
                    self.transform[0],
                    self.transform[1],
                    self.transform[2],
                    -self.transform[3],
                    -self.transform[4],
                    -self.transform[5],
                )
                self.trans_inv = ~self.transform  # pylint: disable=invalid-unary-operand-type

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
                    logging.debug("image contains %d nodata values ", np.sum(self.mask == 0))

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
