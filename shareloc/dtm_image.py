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
Image DTM Image to handle DTM image data.
Inherits from Image Class with DTM particularities.
"""

# Standard imports
import logging

# Third party imports
import numpy as np
from rasterio.fill import fillnodata

# Shareloc imports
from shareloc.image import Image


class DTMImage(Image):
    """
    class DTM  Image to handle DTM image data
    Inherits from Image Class with DTM particularities
    """

    def __init__(
        self,
        image_path,
        read_data=False,
        datum=None,
        roi=None,
        roi_is_in_physical_space=False,
        fill_nodata="rio_fillnodata",
        fill_value=None,
    ):
        """
        constructor
        :param image_path : image path
        :type image_path  : string or None
        :param read_data  : read image data
        :type read_data  : bool
        :param datum  :  datum "geoid" or "ellipsoid", datum is auto identified from babel header if image format is BSQ
           otherwise if None datum is set to "geoid"
        :type datum  : str
        :param roi  : region of interest [row_min,col_min,row_max,col_max] or [xmin,y_min,x_max,y_max] if
             roi_is_in_physical_space activated
        :type roi  : list
        :param roi_is_in_physical_space  : roi value in physical space
        :type roi_is_in_physical_space  : bool
        :param fill_nodata  fill_nodata strategy in None/'constant'/'min'/'median'/'max'/'mean'/'rio_fillnodata'/
        :type fill_nodata  : str
        :param fill_value  fill value for constant strategy. fill value is used for 'roi_fillnodata' residuals nodata,
        if None 'min' is used
        :type fill_value  : float
        """

        super().__init__(image_path, read_data=read_data, roi=roi, roi_is_in_physical_space=roi_is_in_physical_space)
        if datum is None:
            self.datum = "geoid"
        else:
            self.datum = datum

        self.stats = {}
        if read_data:
            if self.mask is not None:
                valid_data = self.data[self.mask[:, :] == 255]
            else:
                valid_data = self.data
            self.stats["min"] = valid_data.min()
            self.stats["max"] = valid_data.max()
            self.stats["mean"] = valid_data.mean()
            self.stats["median"] = np.median(valid_data)
        if fill_nodata is not None:
            self.fill_nodata(strategy=fill_nodata, fill_value=fill_value)

    def fill_nodata(self, strategy="rio_fillnodata", max_search_distance=100.0, smoothing_iterations=0, fill_value=0.0):
        """
        fill nodata in DTM image

        :param strategy: fill strategy ('constant'/'min'/'median'/'max'/'mean'/'rio_fillnodata'/)
        :type strategy: str
        :param max_search_distance: fill max_search_distance
        :type max_search_distance: float
        :param smoothing_iterations: smoothing_iterations
        :type smoothing_iterations: int
        :param fill_value  fill value for constant strategy. fill value is used for 'roi_fillnodata' residuals nodata,
        if None 'min' is used
        :type fill_value  : float
        """
        if self.mask is not None:
            if strategy in self.stats:
                self.data[self.mask[:, :] == 0] = self.stats[strategy]
            elif strategy == "rio_fillnodata":
                self.data = fillnodata(self.data, self.mask[:, :], max_search_distance, smoothing_iterations)
                if np.sum(self.data[self.mask[:, :] == 0] == self.nodata) != 0:
                    if fill_value is None:
                        fill_value = self.stats["min"]
                    logging.warning("not all nodata have been filled, fill with %d", fill_value)
                    self.data[self.data[:, :] == self.nodata] = fill_value
            elif strategy == "constant":
                self.data[self.mask[:, :] == 0] = fill_value
            else:
                logging.warning("fill nodata strategy not available")
        else:
            logging.debug("no nodata mask has been defined")
