#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).
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
This module contains the RPC formats handling to instantiate RPC models.
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
from rasterio.errors import RasterioIOError


def rpc_reader(geomodel_path: str, topleftconvention: bool = True) -> Dict:
    """
    Load from any RPC (auto identify driver)
    from filename (dimap, ossim kwl, geotiff)

    TODO: topleftconvention always to True, set a standard and remove the option

    topleftconvention boolean: [0,0] position
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)

    :param geomodel_path geomodel filename
    :return rpc dict filled with parameters
    """

    rpc_params = rpc_reader_via_rasterio(geomodel_path, topleftconvention)
    if rpc_params is not None:
        return rpc_params

    # If ends with XML --> DIMAP
    if basename(geomodel_path.upper()).endswith("XML"):
        dimap_version = identify_dimap(geomodel_path)
        if dimap_version is not None:
            if float(dimap_version) < 2.0:
                return rpc_reader_dimap_v1(geomodel_path, topleftconvention)
            if float(dimap_version) >= 2.0:
                return rpc_reader_dimap_v23(geomodel_path, topleftconvention)
    ossim_model = identify_ossim_kwl(geomodel_path)
    if ossim_model is not None:
        return rpc_reader_ossim_kwl(geomodel_path, topleftconvention)

    raise ValueError("can not read rpc file")  # noqa: B904


def rpc_reader_dimap_v23(geomodel_path: str, topleftconvention: bool = True) -> Dict:
    """
    Load from Dimap v2 and V3

    :param geomodel_path: path to geomodel
    :param topleftconvention: [0,0] position
    :type topleftconvention: boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)

    :returns rpc_dict
    """

    if not basename(geomodel_path).upper().endswith("XML"):
        raise ValueError("dimap must ends with .xml")

    xmldoc = minidom.parse(geomodel_path)

    mtd = xmldoc.getElementsByTagName("Metadata_Identification")
    version = mtd[0].getElementsByTagName("METADATA_FORMAT")[0].attributes.items()[0][1]
    rpc_params = {"driver_type": "dimap_v" + version}
    global_rfm = xmldoc.getElementsByTagName("Global_RFM")[0]
    normalisation_coeffs = global_rfm.getElementsByTagName("RFM_Validity")[0]
    rpc_params["offset_row"] = float(normalisation_coeffs.getElementsByTagName("LINE_OFF")[0].firstChild.data)
    rpc_params["offset_col"] = float(normalisation_coeffs.getElementsByTagName("SAMP_OFF")[0].firstChild.data)
    if float(version) >= 3:
        direct_coeffs = global_rfm.getElementsByTagName("ImagetoGround_Values")[0]
        inverse_coeffs = global_rfm.getElementsByTagName("GroundtoImage_Values")[0]

        rpc_params["num_x"] = [
            float(direct_coeffs.getElementsByTagName(f"LON_NUM_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_x"] = [
            float(direct_coeffs.getElementsByTagName(f"LON_DEN_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["num_y"] = [
            float(direct_coeffs.getElementsByTagName(f"LAT_NUM_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_y"] = [
            float(direct_coeffs.getElementsByTagName(f"LAT_DEN_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["offset_col"] -= 0.5
        rpc_params["offset_row"] -= 0.5

    else:
        direct_coeffs = global_rfm.getElementsByTagName("Direct_Model")[0]
        inverse_coeffs = global_rfm.getElementsByTagName("Inverse_Model")[0]

        rpc_params["num_x"] = [
            float(direct_coeffs.getElementsByTagName(f"SAMP_NUM_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_x"] = [
            float(direct_coeffs.getElementsByTagName(f"SAMP_DEN_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["num_y"] = [
            float(direct_coeffs.getElementsByTagName(f"LINE_NUM_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_y"] = [
            float(direct_coeffs.getElementsByTagName(f"LINE_DEN_COEFF_{index}")[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["offset_col"] -= 1.0
        rpc_params["offset_row"] -= 1.0

    rpc_params["num_col"] = [
        float(inverse_coeffs.getElementsByTagName(f"SAMP_NUM_COEFF_{index}")[0].firstChild.data)
        for index in range(1, 21)
    ]
    rpc_params["den_col"] = [
        float(inverse_coeffs.getElementsByTagName(f"SAMP_DEN_COEFF_{index}")[0].firstChild.data)
        for index in range(1, 21)
    ]
    rpc_params["num_row"] = [
        float(inverse_coeffs.getElementsByTagName(f"LINE_NUM_COEFF_{index}")[0].firstChild.data)
        for index in range(1, 21)
    ]
    rpc_params["den_row"] = [
        float(inverse_coeffs.getElementsByTagName(f"LINE_DEN_COEFF_{index}")[0].firstChild.data)
        for index in range(1, 21)
    ]

    rpc_params["scale_col"] = float(normalisation_coeffs.getElementsByTagName("SAMP_SCALE")[0].firstChild.data)
    rpc_params["scale_row"] = float(normalisation_coeffs.getElementsByTagName("LINE_SCALE")[0].firstChild.data)
    rpc_params["offset_alt"] = float(normalisation_coeffs.getElementsByTagName("HEIGHT_OFF")[0].firstChild.data)
    rpc_params["scale_alt"] = float(normalisation_coeffs.getElementsByTagName("HEIGHT_SCALE")[0].firstChild.data)
    rpc_params["offset_x"] = float(normalisation_coeffs.getElementsByTagName("LONG_OFF")[0].firstChild.data)
    rpc_params["scale_x"] = float(normalisation_coeffs.getElementsByTagName("LONG_SCALE")[0].firstChild.data)
    rpc_params["offset_y"] = float(normalisation_coeffs.getElementsByTagName("LAT_OFF")[0].firstChild.data)
    rpc_params["scale_y"] = float(normalisation_coeffs.getElementsByTagName("LAT_SCALE")[0].firstChild.data)
    # If top left convention, 0.5 pixel shift added on col/row offsets
    if topleftconvention:
        rpc_params["offset_col"] += 0.5
        rpc_params["offset_row"] += 0.5
    return rpc_params


def rpc_reader_dimap_v1(geomodel_path: str, topleftconvention: bool = True) -> Dict:
    """
    Load from dimap v1

    ** Deprecated, to clean ? **

    :param geomodel_path: path to geomodel
    :param topleftconvention: [0,0] position
    :type topleftconvention: boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
    """

    if not basename(geomodel_path).upper().endswith("XML"):
        raise ValueError("dimap must ends with .xml")

    xmldoc = minidom.parse(geomodel_path)

    mtd = xmldoc.getElementsByTagName("Metadata_Identification")
    version = mtd[0].getElementsByTagName("METADATA_PROFILE")[0].attributes.items()[0][1]
    rpc_params = {"driver_type": "dimap_v" + version}

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

    rpc_params["offset_col"] = offset_col
    rpc_params["scale_col"] = scale_col
    rpc_params["offset_row"] = offset_row
    rpc_params["scale_row"] = scale_row
    rpc_params["offset_alt"] = offset_alt
    rpc_params["scale_alt"] = scale_alt
    rpc_params["offset_x"] = offset_lon
    rpc_params["scale_x"] = scale_lon
    rpc_params["offset_y"] = offset_lat
    rpc_params["scale_y"] = scale_lat
    rpc_params["num_x"] = coeff_lon[0:20]
    rpc_params["den_x"] = coeff_lon[20::]
    rpc_params["num_y"] = coeff_lat[0:20]
    rpc_params["den_y"] = coeff_lat[20::]
    rpc_params["num_col"] = coeff_col[0:20]
    rpc_params["den_col"] = coeff_col[20::]
    rpc_params["num_row"] = coeff_lig[0:20]
    rpc_params["den_row"] = coeff_lig[20::]
    rpc_params["offset_col"] -= 1.0
    rpc_params["offset_row"] -= 1.0
    # If top left convention, 0.5 pixel shift added on col/row offsets
    if topleftconvention:
        rpc_params["offset_col"] += 0.5
        rpc_params["offset_row"] += 0.5
    return rpc_params


def rpc_reader_ossim_kwl(geomodel_path: str, topleftconvention: bool = True) -> Dict:
    """
    Load from a geom file

    :param geomodel_path: path to geomodel
    :param topleftconvention: [0,0] position
    :type topleftconvention: boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
    """
    # OSSIM keyword list
    rpc_params = {"driver_type": "ossim_kwl"}
    with open(geomodel_path, "r", encoding="utf-8") as ossim_file:
        content = ossim_file.readlines()

    geom_dict = {}
    for line in content:
        (key, val) = line.split(": ")
        geom_dict[key] = val.rstrip()

    rpc_params["den_row"] = [np.nan] * 20
    rpc_params["num_row"] = [np.nan] * 20
    rpc_params["den_col"] = [np.nan] * 20
    rpc_params["num_col"] = [np.nan] * 20
    for index in range(0, 20):
        axis = "line"
        num_den = "den"
        key = f"{axis}_{num_den}_coeff_{index:02d}"
        rpc_params["den_row"][index] = float(geom_dict[key])
        num_den = "num"
        key = f"{axis}_{num_den}_coeff_{index:02d}"
        rpc_params["num_row"][index] = float(geom_dict[key])
        axis = "samp"
        key = f"{axis}_{num_den}_coeff_{index:02d}"
        rpc_params["num_col"][index] = float(geom_dict[key])
        num_den = "den"
        key = f"{axis}_{num_den}_coeff_{index:02d}"
        rpc_params["den_col"][index] = float(geom_dict[key])
    rpc_params["offset_col"] = float(geom_dict["samp_off"])
    rpc_params["scale_col"] = float(geom_dict["samp_scale"])
    rpc_params["offset_row"] = float(geom_dict["line_off"])
    rpc_params["scale_row"] = float(geom_dict["line_scale"])
    rpc_params["offset_alt"] = float(geom_dict["height_off"])
    rpc_params["scale_alt"] = float(geom_dict["height_scale"])
    rpc_params["offset_x"] = float(geom_dict["long_off"])
    rpc_params["scale_x"] = float(geom_dict["long_scale"])
    rpc_params["offset_y"] = float(geom_dict["lat_off"])
    rpc_params["scale_y"] = float(geom_dict["lat_scale"])
    # inverse coeff are not defined
    rpc_params["num_x"] = None
    rpc_params["den_x"] = None
    rpc_params["num_y"] = None
    rpc_params["den_y"] = None
    # If top left convention, 0.5 pixel shift added on col/row offsets
    if topleftconvention:
        rpc_params["offset_col"] += 0.5
        rpc_params["offset_row"] += 0.5
    return rpc_params


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


def parse_coeff_line(coeff_str):
    """
    split str coef to float list

    :param coeff_str: line coef
    :type coeff_str: str
    :return: coeff list
    :rtype: list()
    """
    return [float(el) for el in coeff_str.split()]


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


def rpc_reader_via_rasterio(geomodel_path, topleftconvention=True) -> Dict:
    """
    Load via rasterio RPC object

    :param geomodel_path: path to geomodel
    :param topleftconvention: [0,0] position
    :type topleftconvention: boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
    """
    try:
        with rio.open(geomodel_path, "r") as src:
            rpcs = src.rpcs  # pas de coef direct
    except RasterioIOError as rio_error:
        logging.debug("%s can not be read by rasterio", geomodel_path)
        logging.debug("     Erreur Rasterio : %s", rio_error)
        return None

    if not rpcs:
        logging.debug("%s does not contains RPCs readable by rasterio ", geomodel_path)
        return None

    rpcs = rpcs.to_dict()

    rpc_params = {}
    rpc_params["offset_alt"] = rpcs["height_off"]
    rpc_params["scale_alt"] = rpcs["height_scale"]
    rpc_params["offset_y"] = rpcs["lat_off"]
    rpc_params["scale_y"] = rpcs["lat_scale"]
    rpc_params["den_row"] = rpcs["line_den_coeff"]
    rpc_params["num_row"] = rpcs["line_num_coeff"]
    rpc_params["offset_row"] = rpcs["line_off"]
    rpc_params["scale_row"] = rpcs["line_scale"]
    rpc_params["offset_x"] = rpcs["long_off"]
    rpc_params["scale_x"] = rpcs["long_scale"]
    rpc_params["den_col"] = rpcs["samp_den_coeff"]
    rpc_params["num_col"] = rpcs["samp_num_coeff"]
    rpc_params["offset_col"] = rpcs["samp_off"]
    rpc_params["scale_col"] = rpcs["samp_scale"]
    rpc_params["num_x"] = None
    rpc_params["den_x"] = None
    rpc_params["num_y"] = None
    rpc_params["den_y"] = None
    rpc_params["driver_type"] = "rasterio_rpc"

    if topleftconvention:
        rpc_params["offset_col"] += 0.5
        rpc_params["offset_row"] += 0.5

    return rpc_params
