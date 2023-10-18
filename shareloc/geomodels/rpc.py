#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2022 Centre National d'Etudes Spatiales (CNES).
# Copyright (c) 2023 CS GROUP - France, https://csgroup.eu
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
This module contains the RPC class corresponding to the RPC models.
RPC models covered are : DIMAP V1, DIMAP V2, DIMAP V3, ossim (geom file), geotiff.
"""

# Standard imports
import logging
from os.path import basename
from typing import Dict
from xml.dom import minidom

# Third party imports
import numpy as np
import rasterio as rio
from numba import config, njit, prange

# Shareloc imports
from shareloc.geomodels.geomodel import GeoModel
from shareloc.geomodels.geomodel_template import GeoModelTemplate
from shareloc.proj_utils import coordinates_conversion

# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"


def parse_coeff_line(coeff_str):
    """
    split str coef to float list

    :param coeff_str: line coef
    :type coeff_str: str
    :return: coeff list
    :rtype: list()
    """
    return [float(el) for el in coeff_str.split()]


def identify_dimap(xml_file):
    """
    parse xml file to identify dimap and its version

    :param xml_file: dimap rpc file
    :type xml_file: str
    :return: dimap info : dimap_version and None if not an dimap file
    :rtype: str
    """
    try:
        xmldoc = minidom.parse(xml_file)
        mtd = xmldoc.getElementsByTagName("Metadata_Identification")
        mtd_format = mtd[0].getElementsByTagName("METADATA_FORMAT")[0].firstChild.data
        if mtd_format == "DIMAP_PHR":
            version_tag = "METADATA_PROFILE"
        else:
            version_tag = "METADATA_FORMAT"
        version = mtd[0].getElementsByTagName(version_tag)[0].attributes.items()[0][1]
    except Exception:  # pylint: disable=broad-except
        return None

    return version


def identify_ossim_kwl(ossim_kwl_file):
    """
    parse geom file to identify if it is an ossim model

    :param ossim_kwl_fil : ossim keyword list file
    :type ossim_kwl_file: str
    :return: ossimmodel or None if not an ossim kwl file
    :rtype: str
    """
    try:
        with open(ossim_kwl_file, encoding="utf-8") as ossim_file:
            content = ossim_file.readlines()
        geom_dict = {}
        for line in content:
            (key, val) = line.split(": ")
            geom_dict[key] = val.rstrip()
        if "type" in geom_dict:
            if geom_dict["type"].strip().startswith("ossim"):
                return geom_dict["type"].strip()
        return None
    except Exception:  # pylint: disable=broad-except
        return None


def identify_geotiff_rpc(image_filename):
    """
    read image file to identify if it is a geotiff which contains RPCs

    :param image_filename: image_filename
    :type image_filename: str
    :return: rpc info, rpc dict or None  if not a geotiff with rpc
    :rtype: str
    """
    try:
        dataset = rio.open(image_filename)
        rpc_dict = dataset.tags(ns="RPC")
        return rpc_dict
    except Exception:  # pylint: disable=broad-except
        return None


# pylint: disable=no-member


@GeoModel.register("RPC")
class RPC(GeoModelTemplate):
    """
    RPC class including direct and inverse localization instance methods
    """

    # gitlab issue #61
    # pylint: disable=too-many-instance-attributes
    def __init__(self, geomodel_path: str):
        # Instanciate GeoModelTemplate generic init with shared parameters
        super().__init__(geomodel_path)
        self.type = "RPC"

        # initiate epsg and datum, can overriden by rpc_params
        self.epsg = None
        self.datum = None

        # RPC parameters as Dict
        self.rpc_params: Dict = {}

        # RPC parameters are load from geomodel_path to rpc params
        self.load()

        # set class parameters from rpc_params (epsg and datum can be overriden)
        for key, value in self.rpc_params.items():
            setattr(self, key, value)

        if self.epsg is None:
            self.epsg = 4326
        if self.datum is None:
            self.datum = "ellipsoid"

        self.lim_extrapol = 1.0001

        # Monomes seems not used in shareloc code: Clean ?
        # Each monome: c[0]*X**c[1]*Y**c[2]*Z**c[3]
        self.monomes = np.array(
            [
                [1, 0, 0, 0],
                [1, 1, 0, 0],
                [1, 0, 1, 0],
                [1, 0, 0, 1],
                [1, 1, 1, 0],
                [1, 1, 0, 1],
                [1, 0, 1, 1],
                [1, 2, 0, 0],
                [1, 0, 2, 0],
                [1, 0, 0, 2],
                [1, 1, 1, 1],
                [1, 3, 0, 0],
                [1, 1, 2, 0],
                [1, 1, 0, 2],
                [1, 2, 1, 0],
                [1, 0, 3, 0],
                [1, 0, 1, 2],
                [1, 2, 0, 1],
                [1, 0, 2, 1],
                [1, 0, 0, 3],
            ]
        )

        # monomial coefficients of 1st variable derivative
        self.monomes_deriv_1 = np.array(
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [1, 0, 1, 0],
                [1, 0, 0, 1],
                [0, 0, 1, 1],
                [2, 1, 0, 0],
                [0, 0, 2, 0],
                [0, 0, 0, 2],
                [1, 0, 1, 1],
                [3, 2, 0, 0],
                [1, 0, 2, 0],
                [1, 0, 0, 2],
                [2, 1, 1, 0],
                [0, 0, 3, 0],
                [0, 0, 1, 2],
                [2, 1, 0, 1],
                [0, 0, 2, 1],
                [0, 0, 0, 3],
            ]
        )

        # monomial coefficients of 2nd variable derivative
        self.monomes_deriv_2 = np.array(
            [
                [0, 0, 0, 0],
                [0, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 0, 1],
                [1, 1, 0, 0],
                [0, 1, 0, 1],
                [1, 0, 0, 1],
                [0, 2, 0, 0],
                [2, 0, 1, 0],
                [0, 0, 0, 2],
                [1, 1, 0, 1],
                [0, 3, 0, 0],
                [2, 1, 1, 0],
                [0, 1, 0, 2],
                [1, 2, 0, 0],
                [3, 0, 2, 0],
                [1, 0, 0, 2],
                [0, 2, 0, 1],
                [2, 0, 1, 1],
                [0, 0, 0, 3],
            ]
        )

        self.inverse_coefficient = False
        self.direct_coefficient = False

        # pylint: disable=access-member-before-definition
        if self.num_col:
            self.inverse_coefficient = True
            self.num_col = np.array(self.num_col)
            self.den_col = np.array(self.den_col)
            self.num_row = np.array(self.num_row)
            self.den_row = np.array(self.den_row)

        # pylint: disable=access-member-before-definition
        if self.num_x:
            self.direct_coefficient = True
            self.num_x = np.array(self.num_x)
            self.den_x = np.array(self.den_x)
            self.num_y = np.array(self.num_y)
            self.den_y = np.array(self.den_y)

        self.alt_minmax = [self.offset_alt - self.scale_alt, self.offset_alt + self.scale_alt]
        self.col0 = self.offset_col - self.scale_col
        self.colmax = self.offset_col + self.scale_col
        self.row0 = self.offset_row - self.scale_row
        self.rowmax = self.offset_row + self.scale_row

    def load(self):
        """
        Load from any RPC (auto identify driver)
        from filename (dimap, ossim kwl, geotiff)

        TODO: topleftconvention always to True, set a standard and remove the option

        topleftconvention boolean: [0,0] position
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        # Set topleftconvention (keeping historic option): to clean
        topleftconvention = True

        # If ends with XML --> DIMAP
        if basename(self.geomodel_path.upper()).endswith("XML"):
            dimap_version = identify_dimap(self.geomodel_path)
            if dimap_version is not None:
                if float(dimap_version) < 2.0:
                    self.load_dimap_v1(topleftconvention)
                if float(dimap_version) >= 2.0:
                    self.load_dimap_coeff(topleftconvention)
        else:
            # If not DIMAP, identify ossim
            ossim_model = identify_ossim_kwl(self.geomodel_path)
            if ossim_model is not None:
                self.load_ossim_kwl(topleftconvention)
            else:
                # Otherwise, check if RPC is in geotif
                geotiff_rpc_dict = identify_geotiff_rpc(self.geomodel_path)
                if geotiff_rpc_dict is not None:
                    self.load_geotiff(topleftconvention)
                else:
                    # Raise error if no file recognized.
                    raise ValueError("can not read rpc file")

    def load_dimap_coeff(self, topleftconvention=True):
        """
        Load from Dimap v2 and V3

        :param topleftconvention: [0,0] position
        :type topleftconvention: boolean
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """

        if not basename(self.geomodel_path).upper().endswith("XML"):
            raise ValueError("dimap must ends with .xml")

        xmldoc = minidom.parse(self.geomodel_path)

        mtd = xmldoc.getElementsByTagName("Metadata_Identification")
        version = mtd[0].getElementsByTagName("METADATA_FORMAT")[0].attributes.items()[0][1]
        self.rpc_params["driver_type"] = "dimap_v" + version
        global_rfm = xmldoc.getElementsByTagName("Global_RFM")[0]
        normalisation_coeffs = global_rfm.getElementsByTagName("RFM_Validity")[0]
        self.rpc_params["offset_row"] = float(normalisation_coeffs.getElementsByTagName("LINE_OFF")[0].firstChild.data)
        self.rpc_params["offset_col"] = float(normalisation_coeffs.getElementsByTagName("SAMP_OFF")[0].firstChild.data)
        if float(version) >= 3:
            direct_coeffs = global_rfm.getElementsByTagName("ImagetoGround_Values")[0]
            inverse_coeffs = global_rfm.getElementsByTagName("GroundtoImage_Values")[0]

            self.rpc_params["num_x"] = [
                float(direct_coeffs.getElementsByTagName(f"LON_NUM_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["den_x"] = [
                float(direct_coeffs.getElementsByTagName(f"LON_DEN_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["num_y"] = [
                float(direct_coeffs.getElementsByTagName(f"LAT_NUM_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["den_y"] = [
                float(direct_coeffs.getElementsByTagName(f"LAT_DEN_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["offset_col"] -= 0.5
            self.rpc_params["offset_row"] -= 0.5

        else:
            direct_coeffs = global_rfm.getElementsByTagName("Direct_Model")[0]
            inverse_coeffs = global_rfm.getElementsByTagName("Inverse_Model")[0]

            self.rpc_params["num_x"] = [
                float(direct_coeffs.getElementsByTagName(f"SAMP_NUM_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["den_x"] = [
                float(direct_coeffs.getElementsByTagName(f"SAMP_DEN_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["num_y"] = [
                float(direct_coeffs.getElementsByTagName(f"LINE_NUM_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["den_y"] = [
                float(direct_coeffs.getElementsByTagName(f"LINE_DEN_COEFF_{index}")[0].firstChild.data)
                for index in range(1, 21)
            ]
            self.rpc_params["offset_col"] -= 1.0
            self.rpc_params["offset_row"] -= 1.0

        self.rpc_params["num_col"] = [
            float(inverse_coeffs.getElementsByTagName(f"SAMP_NUM_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        self.rpc_params["den_col"] = [
            float(inverse_coeffs.getElementsByTagName(f"SAMP_DEN_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        self.rpc_params["num_row"] = [
            float(inverse_coeffs.getElementsByTagName(f"LINE_NUM_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        self.rpc_params["den_row"] = [
            float(inverse_coeffs.getElementsByTagName(f"LINE_DEN_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]

        self.rpc_params["scale_col"] = float(normalisation_coeffs.getElementsByTagName("SAMP_SCALE")[0].firstChild.data)
        self.rpc_params["scale_row"] = float(normalisation_coeffs.getElementsByTagName("LINE_SCALE")[0].firstChild.data)
        self.rpc_params["offset_alt"] = float(
            normalisation_coeffs.getElementsByTagName("HEIGHT_OFF")[0].firstChild.data
        )
        self.rpc_params["scale_alt"] = float(
            normalisation_coeffs.getElementsByTagName("HEIGHT_SCALE")[0].firstChild.data
        )
        self.rpc_params["offset_x"] = float(normalisation_coeffs.getElementsByTagName("LONG_OFF")[0].firstChild.data)
        self.rpc_params["scale_x"] = float(normalisation_coeffs.getElementsByTagName("LONG_SCALE")[0].firstChild.data)
        self.rpc_params["offset_y"] = float(normalisation_coeffs.getElementsByTagName("LAT_OFF")[0].firstChild.data)
        self.rpc_params["scale_y"] = float(normalisation_coeffs.getElementsByTagName("LAT_SCALE")[0].firstChild.data)
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            self.rpc_params["offset_col"] += 0.5
            self.rpc_params["offset_row"] += 0.5

    def load_dimap_v1(self, topleftconvention=True):
        """
        Load from dimap v1

        ** Deprecated, to clean ? **

        :param topleftconvention: [0,0] position
        :type topleftconvention: boolean
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """

        if not basename(self.geomodel_path).upper().endswith("XML"):
            raise ValueError("dimap must ends with .xml")

        xmldoc = minidom.parse(self.geomodel_path)

        mtd = xmldoc.getElementsByTagName("Metadata_Identification")
        version = mtd[0].getElementsByTagName("METADATA_PROFILE")[0].attributes.items()[0][1]
        self.rpc_params = {"driver_type": "dimap_v" + version}

        global_rfm = xmldoc.getElementsByTagName("Global_RFM")
        rfm_validity = xmldoc.getElementsByTagName("RFM_Validity")
        coeff_lon = [float(el) for el in global_rfm[0].getElementsByTagName("F_LON")[0].firstChild.data.split()]
        coeff_lat = [float(el) for el in global_rfm[0].getElementsByTagName("F_LAT")[0].firstChild.data.split()]
        coeff_col = [float(el) for el in global_rfm[0].getElementsByTagName("F_COL")[0].firstChild.data.split()]
        coeff_lig = [float(el) for el in global_rfm[0].getElementsByTagName("F_ROW")[0].firstChild.data.split()]

        scale_lon = float(rfm_validity[0].getElementsByTagName("Lon")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_lon = float(rfm_validity[0].getElementsByTagName("Lon")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_lat = float(rfm_validity[0].getElementsByTagName("Lat")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_lat = float(rfm_validity[0].getElementsByTagName("Lat")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_alt = float(rfm_validity[0].getElementsByTagName("Alt")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_alt = float(rfm_validity[0].getElementsByTagName("Alt")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_col = float(rfm_validity[0].getElementsByTagName("Col")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_col = float(rfm_validity[0].getElementsByTagName("Col")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_row = float(rfm_validity[0].getElementsByTagName("Row")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_row = float(rfm_validity[0].getElementsByTagName("Row")[0].getElementsByTagName("B")[0].firstChild.data)

        self.rpc_params["offset_col"] = offset_col
        self.rpc_params["scale_col"] = scale_col
        self.rpc_params["offset_row"] = offset_row
        self.rpc_params["scale_row"] = scale_row
        self.rpc_params["offset_alt"] = offset_alt
        self.rpc_params["scale_alt"] = scale_alt
        self.rpc_params["offset_x"] = offset_lon
        self.rpc_params["scale_x"] = scale_lon
        self.rpc_params["offset_y"] = offset_lat
        self.rpc_params["scale_y"] = scale_lat
        self.rpc_params["num_x"] = coeff_lon[0:20]
        self.rpc_params["den_x"] = coeff_lon[20::]
        self.rpc_params["num_y"] = coeff_lat[0:20]
        self.rpc_params["den_y"] = coeff_lat[20::]
        self.rpc_params["num_col"] = coeff_col[0:20]
        self.rpc_params["den_col"] = coeff_col[20::]
        self.rpc_params["num_row"] = coeff_lig[0:20]
        self.rpc_params["den_row"] = coeff_lig[20::]
        self.rpc_params["offset_col"] -= 1.0
        self.rpc_params["offset_row"] -= 1.0
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            self.rpc_params["offset_col"] += 0.5
            self.rpc_params["offset_row"] += 0.5

    def load_geotiff(self, topleftconvention=True):
        """
        Load from a geotiff image file


        :param topleftconvention: [0,0] position
        :type topleftconvention: boolean
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        dataset = rio.open(self.geomodel_path)
        rpc_dict = dataset.tags(ns="RPC")
        if not rpc_dict:
            logging.error("%s does not contains RPCs ", self.geomodel_path)
            raise ValueError
        self.rpc_params = {
            "den_row": parse_coeff_line(rpc_dict["LINE_DEN_COEFF"]),
            "num_row": parse_coeff_line(rpc_dict["LINE_NUM_COEFF"]),
            "num_col": parse_coeff_line(rpc_dict["SAMP_NUM_COEFF"]),
            "den_col": parse_coeff_line(rpc_dict["SAMP_DEN_COEFF"]),
            "offset_col": float(rpc_dict["SAMP_OFF"]),
            "scale_col": float(rpc_dict["SAMP_SCALE"]),
            "offset_row": float(rpc_dict["LINE_OFF"]),
            "scale_row": float(rpc_dict["LINE_SCALE"]),
            "offset_alt": float(rpc_dict["HEIGHT_OFF"]),
            "scale_alt": float(rpc_dict["HEIGHT_SCALE"]),
            "offset_x": float(rpc_dict["LONG_OFF"]),
            "scale_x": float(rpc_dict["LONG_SCALE"]),
            "offset_y": float(rpc_dict["LAT_OFF"]),
            "scale_y": float(rpc_dict["LAT_SCALE"]),
            "num_x": None,
            "den_x": None,
            "num_y": None,
            "den_y": None,
        }
        # inverse coeff are not defined
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            self.rpc_params["offset_col"] += 0.5
            self.rpc_params["offset_row"] += 0.5

    def load_ossim_kwl(self, topleftconvention=True):
        """
        Load from a geom file

        :param topleftconvention: [0,0] position
        :type topleftconvention: boolean
            If False : [0,0] is at the center of the Top Left pixel
            If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        # OSSIM keyword list
        self.rpc_params["driver_type"] = "ossim_kwl"
        with open(self.geomodel_path, "r", encoding="utf-8") as ossim_file:
            content = ossim_file.readlines()

        geom_dict = {}
        for line in content:
            (key, val) = line.split(": ")
            geom_dict[key] = val.rstrip()

        self.rpc_params["den_row"] = [np.nan] * 20
        self.rpc_params["num_row"] = [np.nan] * 20
        self.rpc_params["den_col"] = [np.nan] * 20
        self.rpc_params["num_col"] = [np.nan] * 20
        for index in range(0, 20):
            axis = "line"
            num_den = "den"
            key = f"{axis}_{num_den}_coeff_{index:02d}"
            self.rpc_params["den_row"][index] = float(geom_dict[key])
            num_den = "num"
            key = f"{axis}_{num_den}_coeff_{index:02d}"
            self.rpc_params["num_row"][index] = float(geom_dict[key])
            axis = "samp"
            key = f"{axis}_{num_den}_coeff_{index:02d}"
            self.rpc_params["num_col"][index] = float(geom_dict[key])
            num_den = "den"
            key = f"{axis}_{num_den}_coeff_{index:02d}"
            self.rpc_params["den_col"][index] = float(geom_dict[key])
        self.rpc_params["offset_col"] = float(geom_dict["samp_off"])
        self.rpc_params["scale_col"] = float(geom_dict["samp_scale"])
        self.rpc_params["offset_row"] = float(geom_dict["line_off"])
        self.rpc_params["scale_row"] = float(geom_dict["line_scale"])
        self.rpc_params["offset_alt"] = float(geom_dict["height_off"])
        self.rpc_params["scale_alt"] = float(geom_dict["height_scale"])
        self.rpc_params["offset_x"] = float(geom_dict["long_off"])
        self.rpc_params["scale_x"] = float(geom_dict["long_scale"])
        self.rpc_params["offset_y"] = float(geom_dict["lat_off"])
        self.rpc_params["scale_y"] = float(geom_dict["lat_scale"])
        # inverse coeff are not defined
        self.rpc_params["num_x"] = None
        self.rpc_params["den_x"] = None
        self.rpc_params["num_y"] = None
        self.rpc_params["den_y"] = None
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            self.rpc_params["offset_col"] += 0.5
            self.rpc_params["offset_row"] += 0.5

    def direct_loc_h(self, row, col, alt, fill_nan=False):
        """
        direct localization at constant altitude

        :param row:  line sensor position
        :type row: float or 1D numpy.ndarray dtype=float64
        :param col:  column sensor position
        :type col: float or 1D numpy.ndarray dtype=float64
        :param alt:  altitude
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan: boolean
        :return: ground position (lon,lat,h)
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """
        if not isinstance(col, (list, np.ndarray)):
            col = np.array([col])
            row = np.array([row])

        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != col.shape[0]:
            alt = np.full(col.shape[0], fill_value=alt[0])

        points = np.zeros((col.size, 3))
        filter_nan, points[:, 0], points[:, 1] = self.filter_coordinates(row, col, fill_nan)
        row = row[filter_nan]
        col = col[filter_nan]

        # Direct localization using direct RPC
        if self.direct_coefficient:
            # ground position
            col_norm = (col - self.offset_col) / self.scale_col
            row_norm = (row - self.offset_row) / self.scale_row
            alt_norm = (alt - self.offset_alt) / self.scale_alt

            if np.sum(abs(col_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! column extrapolation in direct localization ")
            if np.sum(abs(row_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! row extrapolation in direct localization ")
            if np.sum(abs(alt_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! alt extrapolation in direct localization ")

            points[filter_nan, 1], points[filter_nan, 0] = compute_rational_function_polynomial(
                col_norm,
                row_norm,
                alt_norm,
                self.num_x,
                self.den_x,
                self.num_y,
                self.den_y,
                self.scale_x,
                self.offset_x,
                self.scale_y,
                self.offset_y,
            )

        # Direct localization using inverse RPC
        else:
            logging.debug("direct localisation from inverse iterative")
            (points[filter_nan, 0], points[filter_nan, 1], points[filter_nan, 2]) = self.direct_loc_inverse_iterative(
                row, col, alt, 10, fill_nan
            )
        points[:, 2] = alt
        return points

    def direct_loc_grid_h(self, row0, col0, steprow, stepcol, nbrow, nbcol, alt):
        """
        calculates a direct loc grid (lat, lon) from the direct RPCs at constant altitude
        TODO: not tested.

        :param row0:  grid origin (row)
        :type row0: int
        :param col0:  grid origin (col)
        :type col0: int
        :param steprow:  grid step (row)
        :type steprow: int
        :param stepcol:  grid step (col)
        :type stepcol: int
        :param nbrow:  grid nb row
        :type nbrow: int
        :param nbcol:  grid nb col
        :type nbcol: int
        :param alt: altitude of the grid
        :type alt: float
        :return: direct localization grid longitude and latitude
        :rtype: Tuple(numpy.array, numpy.array)
        """
        gri_lon = np.zeros((nbrow, nbcol))
        gri_lat = np.zeros((nbrow, nbcol))
        for column in range(int(nbcol)):
            col = col0 + stepcol * column
            for line in range(int(nbrow)):
                row = row0 + steprow * line
                (gri_lon[line, column], gri_lat[line, column], __) = self.direct_loc_h(row, col, alt)
        return (gri_lon, gri_lat)

    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm

        :param row:  line sensor position
        :type row: float
        :param col:  column sensor position
        :type col: float
        :param dtm: dtm intersection model
        :type dtm: shareloc.geofunctions.dtm_intersection
        :return: ground position (lon,lat,h) in dtm coordinates system
        :rtype: numpy.ndarray 2D dimension with (N,3) shape, where N is number of input coordinates
        """
        if isinstance(col, (list, np.ndarray)):
            points_nb = len(col)
        else:
            points_nb = 1
            row = np.array([row])
            col = np.array([col])
        direct_dtm = np.zeros((points_nb, 3))

        diff_alti_min, diff_alti_max = dtm.get_alt_offset(self.epsg)
        # print("min {} max {}".format(dtm.Zmin,dtm.Zmax))
        (min_dtm, max_dtm) = (dtm.alt_min - 1.0 + diff_alti_min, dtm.alt_max + 1.0 + diff_alti_max)
        if min_dtm < self.offset_alt - self.scale_alt:
            logging.debug("minimum dtm value is outside RPC validity domain, extrapolation will be done")
        if max_dtm > self.offset_alt + self.scale_alt:
            logging.debug("maximum dtm value is outside RPC validity domain, extrapolation will be done")
        los = self.los_extrema(row, col, min_dtm, max_dtm, epsg=dtm.epsg)
        for i in range(points_nb):
            los_i = los[2 * i : 2 * i + 2, :]
            (__, __, position_cube, alti, los_index) = dtm.intersect_dtm_cube(los_i)
            if position_cube is not None:
                (__, __, position) = dtm.intersection(los_index, position_cube, alti)
                direct_dtm[i, :] = position
            else:
                position = np.full(3, fill_value=np.nan)
            direct_dtm[i, :] = position
        return direct_dtm

    def inverse_loc(self, lon, lat, alt):
        """
        Inverse localization

        :param lon: longitude position
        :type lon: float or 1D numpy.ndarray dtype=float64
        :param lat: latitude position
        :type lat: float or 1D numpy.ndarray dtype=float64
        :param alt: altitude
        :type alt: float
        :return: sensor position (row, col, alt)
        :rtype: tuple(1D np.array row position, 1D np.array col position, 1D np.array alt)
        """
        if self.inverse_coefficient:
            if not isinstance(lon, (list, np.ndarray)):
                lon = np.array([lon])
                lat = np.array([lat])
            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != lon.shape[0]:
                alt = np.full(lon.shape[0], fill_value=alt[0])

            lon_norm = (lon - self.offset_x) / self.scale_x
            lat_norm = (lat - self.offset_y) / self.scale_y
            alt_norm = (alt - self.offset_alt) / self.scale_alt

            if np.sum(abs(lon_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! longitude extrapolation in inverse localization ")
            if np.sum(abs(lat_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! row extrapolation in inverse localization ")
            if np.sum(abs(alt_norm) > self.lim_extrapol) > 0:
                logging.debug("!!!!! alt extrapolation in inverse localization ")

            row_out, col_out = compute_rational_function_polynomial(
                lon_norm,
                lat_norm,
                alt_norm,
                self.num_col,
                self.den_col,
                self.num_row,
                self.den_row,
                self.scale_col,
                self.offset_col,
                self.scale_row,
                self.offset_row,
            )
        else:
            logging.warning("inverse localisation can't be performed, inverse coefficients have not been defined")
            (col_out, row_out) = (None, None)
        return row_out, col_out, alt

    def filter_coordinates(self, first_coord, second_coord, fill_nan=False, direction="direct"):
        """
        Filter nan input values

        :param first_coord:  first coordinate
        :type first_coord: 1D numpy.ndarray dtype=float64
        :param second_coord:  second coordinate
        :type second_coord: 1D numpy.ndarray dtype=float64
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan: boolean
        :param direction:  direct or inverse localisation
        :type direction: str in ('direct', 'inverse')
        :return: filtered coordinates
        :rtype: list of numpy.array (index of nan, first filtered, second filtered)
        """
        filter_nan = np.logical_not(np.logical_or(np.isnan(first_coord), np.isnan(second_coord)))

        if fill_nan:
            if direction == "direct":
                out_x_nan_value = self.offset_x
                out_y_nan_value = self.offset_y
            else:
                out_x_nan_value = self.offset_col
                out_y_nan_value = self.offset_row
        else:
            out_x_nan_value = np.nan
            out_y_nan_value = np.nan

        x_out = np.full(len(second_coord), out_x_nan_value)
        y_out = np.full(len(second_coord), out_y_nan_value)

        return filter_nan, x_out, y_out

    def compute_loc_inverse_derivates(self, lon, lat, alt):
        """
        Inverse loc partial derivatives analytical compute

        :param lon: longitude coordinate
        :param lat: latitude coordinate
        :param alt: altitude coordinate
        :return: partials derivatives of inverse localization
        :rtype: Tuple(dcol_dlon np.array, dcol_dlat np.array, drow_dlon np.array, drow_dlat np.array)
        """
        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != lon.shape[0]:
            alt = np.full(lon.shape[0], fill_value=alt[0])

        lon_norm = (lon - self.offset_x) / self.scale_x
        lat_norm = (lat - self.offset_y) / self.scale_y
        alt_norm = (alt - self.offset_alt) / self.scale_alt

        dcol_dlon, dcol_dlat, drow_dlon, drow_dlat = compute_loc_inverse_derivates_numba(
            lon_norm,
            lat_norm,
            alt_norm,
            self.num_col,
            self.den_col,
            self.num_row,
            self.den_row,
            self.scale_col,
            self.scale_x,
            self.scale_row,
            self.scale_y,
        )

        return (dcol_dlon, dcol_dlat, drow_dlon, drow_dlat)

    def direct_loc_inverse_iterative(self, row, col, alt, nb_iter_max=10, fill_nan=False):
        """
        Iterative direct localization using inverse RPC

        :param row:  line sensor position
        :type row: float or 1D numpy.ndarray dtype=float64
        :param col:  column sensor position
        :type col: float or 1D numpy.ndarray dtype=float64
        :param alt:  altitude
        :type alt: float
        :param nb_iter_max: max number of iteration
        :type alt: int
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan: boolean
        :return: ground position (lon,lat,h)
        :rtype: list of numpy.array
        """
        if self.inverse_coefficient:
            if not isinstance(row, (list, np.ndarray)):
                col = np.array([col])
                row = np.array([row])

            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != col.shape[0]:
                alt = np.full(col.shape[0], fill_value=alt[0])

            filter_nan, long_out, lat_out = self.filter_coordinates(row, col, fill_nan)
            row = row[filter_nan]
            col = col[filter_nan]
            alt = alt[filter_nan]

            # if all coord contains Nan then return
            if not np.any(filter_nan):
                return long_out, lat_out, alt

            # inverse localization starting from the center of the scene
            lon = np.array([self.offset_x])
            lat = np.array([self.offset_y])
            (row_start, col_start, __) = self.inverse_loc(lon, lat, alt)

            # desired precision in pixels
            eps = 1e-6

            iteration = 0
            # computing the residue between the sensor positions and those estimated by the inverse localization
            delta_col = col - col_start
            delta_row = row - row_start

            # ground coordinates (latitude and longitude) of each point
            lon = np.repeat(lon, delta_col.size)
            lat = np.repeat(lat, delta_col.size)
            # while the required precision is not achieved
            while (np.max(abs(delta_col)) > eps or np.max(abs(delta_row)) > eps) and iteration < nb_iter_max:
                # list of points that require another iteration
                iter_ = np.where((abs(delta_col) > eps) | (abs(delta_row) > eps))[0]

                # partial derivatives
                (dcol_dlon, dcol_dlat, drow_dlon, drow_dlat) = self.compute_loc_inverse_derivates(
                    lon[iter_], lat[iter_], alt[iter_]
                )
                det = dcol_dlon * drow_dlat - drow_dlon * dcol_dlat

                delta_lon = (drow_dlat * delta_col[iter_] - dcol_dlat * delta_row[iter_]) / det
                delta_lat = (-drow_dlon * delta_col[iter_] + dcol_dlon * delta_row[iter_]) / det

                # update ground coordinates
                lon[iter_] += delta_lon
                lat[iter_] += delta_lat

                # inverse localization
                (row_estim, col_estim, __) = self.inverse_loc(lon[iter_], lat[iter_], alt[iter_])

                # updating the residue between the sensor positions and those estimated by the inverse localization
                delta_col[iter_] = col[iter_] - col_estim
                delta_row[iter_] = row[iter_] - row_estim
                iteration += 1

            long_out[filter_nan] = lon
            lat_out[filter_nan] = lat

        else:
            logging.warning("inverse localisation can't be performed, inverse coefficients have not been defined")
            (long_out, lat_out) = (None, None)

        return long_out, lat_out, alt

    def get_alt_min_max(self):
        """
        returns altitudes min and max layers

        :return: alt_min,lat_max
        :rtype: list
        """
        return [self.offset_alt - self.scale_alt / 2.0, self.offset_alt + self.scale_alt / 2.0]

    def los_extrema(self, row, col, alt_min=None, alt_max=None, fill_nan=False, epsg=None):
        """
        compute los extrema

        :param row:  line sensor position
        :type row: float
        :param col:  column sensor position
        :type col: float
        :param alt_min: los alt min
        :type alt_min: float
        :param alt_max: los alt max
        :type alt_max: float
        :param epsg: epsg code of the dtm
        :type epsg: int
        :return: los extrema
        :rtype: numpy.array (2x3)
        """
        extrapolate = False
        if alt_min is None or alt_max is None:
            [los_alt_min, los_alt_max] = self.get_alt_min_max()
        elif alt_min >= self.alt_minmax[0] and alt_max <= self.alt_minmax[1]:
            los_alt_min = alt_min
            los_alt_max = alt_max
        else:
            extrapolate = True
            [los_alt_min, los_alt_max] = self.get_alt_min_max()

        #
        if isinstance(row, (np.ndarray)):
            los_nb = row.shape[0]
            row_array = np.full([los_nb * 2], fill_value=0.0)
            col_array = np.full([los_nb * 2], fill_value=0.0)
            alt_array = np.full([los_nb * 2], fill_value=0.0)
            row_array[0::2] = row
            row_array[1::2] = row
            col_array[0::2] = col
            col_array[1::2] = col
            alt_array[0::2] = los_alt_max
            alt_array[1::2] = los_alt_min
        else:
            los_nb = 1
            row_array = np.array([row, row])
            col_array = np.array([col, col])
            alt_array = np.array([los_alt_max, los_alt_min])
        los_edges = self.direct_loc_h(row_array, col_array, alt_array, fill_nan)
        if extrapolate:
            diff = los_edges[0::2, :] - los_edges[1::2, :]
            delta_alt = diff[:, 2]
            coeff_alt_max = (alt_max - los_edges[1::2, 2]) / delta_alt
            coeff_alt_max = np.tile(coeff_alt_max[:, np.newaxis], (1, 3))
            coeff_alt_min = (alt_min - los_edges[1::2, 2]) / delta_alt
            coeff_alt_min = np.tile(coeff_alt_min[:, np.newaxis], (1, 3))
            los_edges[0::2, :] = los_edges[1::2, :] + diff * coeff_alt_max
            los_edges[1::2, :] = los_edges[1::2, :] + diff * coeff_alt_min
        if epsg is not None and epsg != self.epsg:
            los_edges = coordinates_conversion(los_edges, self.epsg, epsg)
        return los_edges


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def polynomial_equation(xnorm, ynorm, znorm, coeff):
    """
    Compute polynomial equation

    :param xnorm: Normalized longitude (for inverse) or column (for direct) position
    :type xnorm: float 64
    :param ynorm: Normalized latitude (for inverse) or line (for direct) position
    :type ynorm: float 64
    :param znorm: Normalized altitude position
    :type znorm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational
    :rtype: float 64
    """
    rational = (
        coeff[0]
        + coeff[1] * xnorm
        + coeff[2] * ynorm
        + coeff[3] * znorm
        + coeff[4] * xnorm * ynorm
        + coeff[5] * xnorm * znorm
        + coeff[6] * ynorm * znorm
        + coeff[7] * xnorm**2
        + coeff[8] * ynorm**2
        + coeff[9] * znorm**2
        + coeff[10] * xnorm * ynorm * znorm
        + coeff[11] * xnorm**3
        + coeff[12] * xnorm * ynorm**2
        + coeff[13] * xnorm * znorm**2
        + coeff[14] * xnorm**2 * ynorm
        + coeff[15] * ynorm**3
        + coeff[16] * ynorm * znorm**2
        + coeff[17] * xnorm**2 * znorm
        + coeff[18] * ynorm**2 * znorm
        + coeff[19] * znorm**3
    )

    return rational


# pylint: disable=too-many-arguments
@njit(
    "Tuple((f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)",
    parallel=True,
    cache=True,
    fastmath=True,
)
def compute_rational_function_polynomial(
    lon_col_norm,
    lat_row_norm,
    alt_norm,
    num_col,
    den_col,
    num_lin,
    den_lin,
    scale_col,
    offset_col,
    scale_lin,
    offset_lin,
):
    """
    Compute rational function polynomial using numba to reduce calculation time on multiple points.
    useful to compute direct and inverse localization using direct or inverse RPC.

    :param lon_col_norm: Normalized longitude (for inverse) or column (for direct) position
    :type lon_col_norm: 1D np.array dtype np.float 64
    :param lat_row_norm: Normalized latitude (for inverse) or line (for direct) position
    :type lat_row_norm: 1D np.array dtype np.float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: 1D np.array dtype np.float 64
    :param num_col: Column numerator coefficients
    :type num_col: 1D np.array dtype np.float 64
    :param den_col: Column denominator coefficients
    :type den_col: 1D np.array dtype np.float 64
    :param num_lin: Line numerator coefficients
    :type num_lin: 1D np.array dtype np.float 64
    :param den_lin: Line denominator coefficients
    :type den_lin: 1D np.array dtype np.float 64
    :param scale_col: Column scale
    :type scale_col: float 64
    :param offset_col: Column offset
    :type offset_col: float 64
    :param scale_lin: Line scale
    :type scale_lin: float 64
    :param offset_lin: Line offset
    :type offset_lin: float 64
    :return: for inverse localization : sensor position (row, col). for direct localization : ground position (lon, lat)
    :rtype: Tuple(np.ndarray, np.ndarray)
    """
    assert lon_col_norm.shape == alt_norm.shape

    col_lat_out = np.zeros((lon_col_norm.shape[0]), dtype=np.float64)
    row_lon_out = np.zeros((lon_col_norm.shape[0]), dtype=np.float64)

    # pylint: disable=not-an-iterable
    for i in prange(lon_col_norm.shape[0]):
        poly_num_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_col)
        poly_den_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_col)
        poly_num_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_lin)
        poly_den_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_lin)
        col_lat_out[i] = poly_num_col / poly_den_col * scale_col + offset_col
        row_lon_out[i] = poly_num_lin / poly_den_lin * scale_lin + offset_lin

    return row_lon_out, col_lat_out


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, coeff):
    """
    Compute latitude derivative polynomial equation

    :param lon_norm: Normalized longitude position
    :type lon_norm: float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    derivate = (
        coeff[2]
        + coeff[4] * lon_norm
        + coeff[6] * alt_norm
        + 2 * coeff[8] * lat_norm
        + coeff[10] * lon_norm * alt_norm
        + 2 * coeff[12] * lon_norm * lat_norm
        + coeff[14] * lon_norm**2
        + 3 * coeff[15] * lat_norm**2
        + coeff[16] * alt_norm**2
        + 2 * coeff[18] * lat_norm * alt_norm
    )

    return derivate


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, coeff):
    """
    Compute longitude derivative polynomial equation

    :param lon_norm: Normalized longitude position
    :type lon_norm: float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    derivate = (
        coeff[1]
        + coeff[4] * lat_norm
        + coeff[5] * alt_norm
        + 2 * coeff[7] * lon_norm
        + coeff[10] * lat_norm * alt_norm
        + 3 * coeff[11] * lon_norm**2
        + coeff[12] * lat_norm**2
        + coeff[13] * alt_norm**2
        + 2 * coeff[14] * lat_norm * lon_norm
        + 2 * coeff[17] * lon_norm * alt_norm
    )

    return derivate


# pylint: disable=too-many-arguments
@njit(
    "Tuple((f8[:], f8[:], f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)",
    parallel=True,
    cache=True,
    fastmath=True,
)
def compute_loc_inverse_derivates_numba(
    lon_norm, lat_norm, alt_norm, num_col, den_col, num_lin, den_lin, scale_col, scale_lon, scale_lin, scale_lat
):
    """
    Analytically compute the partials derivatives of inverse localization using numba to reduce calculation time on
    multiple points

    :param lon_norm: Normalized longitude position
    :type lon_norm: 1D np.array dtype np.float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: 1D np.array dtype np.float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: 1D np.array dtype np.float 64
    :param num_col: Column numerator coefficients
    :type num_col: 1D np.array dtype np.float 64
    :param den_col: Column denominator coefficients
    :type den_col: 1D np.array dtype np.float 64
    :param num_lin: Line numerator coefficients
    :type num_lin: 1D np.array dtype np.float 64
    :param den_lin: Line denominator coefficients
    :type den_lin: 1D np.array dtype np.float 64
    :param scale_col: Column scale
    :type scale_col: float 64
    :param scale_lon: Geodetic longitude scale
    :type scale_lon: float 64
    :param scale_lin: Line scale
    :type scale_lin: float 64
    :param scale_lat: Geodetic latitude scale
    :type scale_lat: float 64
    :return: partials derivatives of inverse localization
    :rtype: Tuples(dcol_dlon np.array, dcol_dlat np.array, drow_dlon np.array, drow_dlat np.array)
    """
    dcol_dlon = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    dcol_dlat = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    drow_dlon = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    drow_dlat = np.zeros((lon_norm.shape[0]), dtype=np.float64)

    # pylint: disable=not-an-iterable
    for i in prange(lon_norm.shape[0]):
        num_dcol = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        num_dcol_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        num_dcol_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        dcol_dlon[i] = scale_col / scale_lon * (num_dcol_dlon * den_dcol - den_dcol_dlon * num_dcol) / den_dcol**2
        dcol_dlat[i] = scale_col / scale_lat * (num_dcol_dlat * den_dcol - den_dcol_dlat * num_dcol) / den_dcol**2
        drow_dlon[i] = scale_lin / scale_lon * (num_drow_dlon * den_drow - den_drow_dlon * num_drow) / den_drow**2
        drow_dlat[i] = scale_lin / scale_lat * (num_drow_dlat * den_drow - den_drow_dlat * num_drow) / den_drow**2

    return dcol_dlon, dcol_dlat, drow_dlon, drow_dlat
