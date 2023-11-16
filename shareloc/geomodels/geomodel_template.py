#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of shareloc
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
This module contains the coregistration class template.
It contains the structure for all coregistration methods in subclasses and
generic coregistration code to avoid duplication.
"""

# Standard imports
from abc import ABCMeta, abstractmethod

# Global variable for optimization mode (functions in C)
# SHARELOC_OPTIM_GEOMODEL = False

# TODO: Override functions depending on optimization or not

# if(SHARELOC_OPTIM_GEOMODEL == True):
#     GeoModelTemplate.direct_loc_dtm = GeoModelTemplate.direct_loc_dtm_optim


# Standard imports
import logging
from os.path import basename
from typing import Dict
from xml.dom import minidom

# Third party imports
import numpy as np
import rasterio as rio
from numba import config, njit, prange

from shareloc.proj_utils import coordinates_conversion









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












class GeoModelTemplate(object):#metaclass=ABCMeta):
    """
    Class for general specification of a geometric model
    declined in rpc.py and grid.py
    """

    @abstractmethod
    def __init__(self, geomodel_path: str):
        
        """
        Return the geomodel object associated with the geomodel_type
        given in the configuration

        :param geomodel_path: path of the geomodel file to instanciate.
        :type geomodel_path: string
        """
        # geomodel filename path
        self.geomodel_path: str = geomodel_path

        # geomodel type. Set by the subclass
        self.type: str

        # geomodel filename path
        self.rpc_params: Dict = {}

        self.load()

    # Define GeoModelTemplate functions interface

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

