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

"""
Image DTM Image to handle DTM image data.
Inherits from Image Class with DTM particularities.
"""

# Standard imports
import logging

# Third party imports
import numpy as np
from rasterio.fill import fillnodata
from scipy import interpolate

# Shareloc imports
from shareloc.image import Image
from shareloc.proj_utils import (
    coordinates_conversion,
    transform_index_to_physical_point,
    transform_physical_point_to_index,
)


# pylint: disable=invalid-name
class dtm_reader(Image):
    """
    class dtm_reader to handle DTM image data
    Inherits from Image Class with DTM particularities
    """

    # pylint: disable=too-many-arguments
    def __init__(
        self,
        dtm_filename,
        geoid_filename=None,
        roi=None,
        roi_is_in_physical_space=False,
        fill_nodata="rio_fillnodata",
        fill_value=None,
    ):
        """
        constructor

        :param dtm_filename: dtm filename
        :type dtm_filename: string
        :param geoid_filename: geoid filename, if None datum is ellispoid
        :type geoid_filename: string
        :param roi: region of interest [row_min,col_min,row_max,col_max] or [xmin,y_min,x_max,y_max] if
            roi_is_in_physical_space activated
        :type roi: list
        :param roi_is_in_physical_space: roi value in physical space
        :type roi_is_in_physical_space: bool
        :param fill_nodata:  fill_nodata strategy in None/'constant'/'min'/'median'/'max'/'mean'/'rio_fillnodata'/
        :type fill_nodata: str
        :param fill_value:  fill value for constant strategy. fill value is used for 'roi_fillnodata' residuals nodata,
        if None 'min' is used
        :type fill_value: float
        """

        super().__init__(dtm_filename, read_data=True, roi=roi, roi_is_in_physical_space=roi_is_in_physical_space)

        self.dtm_filename = dtm_filename
        self.geoid_filename = geoid_filename

        self.stats = {}

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
        self.alt_data = self.data[:, :].astype("float64")

        if geoid_filename is not None:
            logging.debug("remove geoid height")
            self.grid_row, self.grid_col = np.mgrid[0 : self.nb_rows : 1, 0 : self.nb_columns : 1]
            lat, lon = transform_index_to_physical_point(self.transform, self.grid_row, self.grid_col)
            positions = np.vstack([lon.flatten(), lat.flatten()]).transpose()
            if self.epsg != 4326:
                positions = coordinates_conversion(positions, self.epsg, 4326)
            geoid_height = interpolate_geoid_height(geoid_filename, positions)
            self.alt_data += geoid_height.reshape(lon.shape)

        else:
            logging.debug("no geoid file is given dtm is assumed to be w.r.t ellipsoid")

        self.trans_inv = self.trans_inv.to_gdal()
        self.transform = self.transform.to_gdal()

    def fill_nodata(self, strategy="rio_fillnodata", max_search_distance=100.0, smoothing_iterations=0, fill_value=0.0):
        """
        fill nodata in DTM image

        :param strategy: fill strategy ('constant'/'min'/'median'/'max'/'mean'/'rio_fillnodata'/)
        :type strategy: str
        :param max_search_distance: fill max_search_distance
        :type max_search_distance: float
        :param smoothing_iterations: smoothing_iterations
        :type smoothing_iterations: int
        :param fill_value: fill value for constant strategy. fill value is used for 'roi_fillnodata' residuals nodata,
            if None 'min' is used
        :type fill_value: float
        """
        if self.mask is not None:
            if strategy in self.stats:
                self.data[self.mask[:, :] == 0] = self.stats[strategy]
            elif strategy == "rio_fillnodata":
                self.data = fillnodata(self.data, self.mask[:, :], max_search_distance, smoothing_iterations)
                if np.sum(self.data[self.mask[:, :] == 0] == self.nodata) != 0:
                    if fill_value is None:
                        fill_value = self.stats["min"]
                    logging.info("Shareloc dtm_reader: not all nodata have been filled, fill with %d", fill_value)
                    self.data[self.data[:, :] == self.nodata] = fill_value
            elif strategy == "constant":
                self.data[self.mask[:, :] == 0] = fill_value
            else:
                logging.warning("Shareloc dtm_reader: fill nodata strategy not available")
        else:
            logging.debug("Shareloc dtm_reader: no nodata mask has been defined")


def interpolate_geoid_height(geoid_filename, positions, interpolation_method="linear"):
    """
    terrain to index conversion
    retrieve geoid height above ellispoid

    :param geoid_filename: geoid_filename
    :type geoid_filename: str
    :param positions: geodetic coordinates
    :type positions: 2D numpy array: (number of points, [long coord, lat coord])
    :param interpolation_method: default is 'linear' (interpn interpolation method)
    :type interpolation_method: str
    :return: geoid height
    :rtype: numpy array (number of points)
    """

    geoid_image = Image(geoid_filename, read_data=True)

    # Check longitude overlap is not present, rounding to handle egm2008 with rounded pixel size
    if geoid_image.nb_columns * geoid_image.pixel_size_col - 360 < 10**-8:
        logging.debug("add one pixel overlap on longitudes")
        geoid_image.nb_columns += 1
        # Check if we can add a column
        geoid_image.data = np.column_stack((geoid_image.data[:, :], geoid_image.data[:, 0]))

    # Prepare grid for interpolation
    row_indexes = np.arange(0, geoid_image.nb_rows, 1)
    col_indexes = np.arange(0, geoid_image.nb_columns, 1)
    points = (row_indexes, col_indexes)

    # add modulo lon/lat
    min_lon = geoid_image.origin_col + geoid_image.pixel_size_col / 2
    max_lon = (
        geoid_image.origin_col + geoid_image.nb_columns * geoid_image.pixel_size_col - geoid_image.pixel_size_col / 2
    )
    positions[:, 0] += ((positions[:, 0] + min_lon) < 0) * 360.0
    positions[:, 0] -= ((positions[:, 0] - max_lon) > 0) * 360.0
    if np.any(np.abs(positions[:, 1]) > 90.0):
        raise RuntimeError("Geoid cannot handle latitudes greater than 90 deg.")
    indexes_geoid = transform_physical_point_to_index(geoid_image.trans_inv, positions[:, 1], positions[:, 0])
    return interpolate.interpn(
        points,
        geoid_image.data[:, :],
        indexes_geoid,
        method=interpolation_method,
        bounds_error=False,
        fill_value=np.nan,
    )
