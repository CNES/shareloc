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

import logging
import numpy as np
from affine import Affine
from shareloc.image.readwrite import read_hdbabel_header, read_bsq_grid
from shareloc.image.image import Image


class DTMImage(Image):
    """ class DTM  Image to handle DTM image data """

    def __init__(self, image_path, read_data=False, datum=None):
        if image_path.split(".")[-1] == "c1":
            logging.debug("bsq babel image")
            if image_path is not None:
                # Image path
                self.image_path = image_path

                babel_dict = read_hdbabel_header(image_path)

                # Pixel size
                self.pixel_size_row = babel_dict["pixel_size_y"]
                self.pixel_size_col = babel_dict["pixel_size_x"]

                # Georeferenced coordinates of the upper-left origin
                # babel origin is pixel centered
                self.origin_row = babel_dict["origin_y"] - self.pixel_size_row / 2.0
                self.origin_col = babel_dict["origin_x"] - self.pixel_size_col / 2.0

                # Image size
                self.nb_rows = babel_dict["row_nb"]
                self.nb_columns = babel_dict["column_nb"]

                self.data_type = babel_dict["data_type"]

                # Geo-transform of type Affine with convention :
                # | pixel size col,   row rotation, origin col |
                # | col rotation  , pixel size row, origin row |
                self.transform = Affine(
                    self.pixel_size_col, 0.0, self.origin_col, 0.0, self.pixel_size_row, self.origin_row
                )

                if "Z-M" in babel_dict["gdlib_code"].split(":")[-1]:
                    self.datum = "ellipsoid"
                elif "H-M" in babel_dict["gdlib_code"].split(":")[-1]:
                    self.datum = "geoid"
                else:
                    logging.warning(
                        "unknown datum {} ellipsoid " "will be used", babel_dict["gdlib_code"].split(":")[-1]
                    )
                    self.datum = "ellipsoid"

                self.data = np.zeros((1, self.nb_rows, self.nb_columns), dtype=self.data_type)
                if read_data:
                    # Data of shape (nb band, nb row, nb col)
                    self.data[0, :, :] = read_bsq_grid(self.image_path, self.nb_rows, self.nb_columns, self.data_type)
        else:
            super().__init__(image_path, read_data=read_data)
            if datum is None:
                self.datum = "geoid"
            else:
                self.datum = datum
