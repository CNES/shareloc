#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2025 Centre National d'Etudes Spatiales (CNES).
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
This module contains the tools to write RPCs
"""

import json

# Standard imports
import logging
import os
from typing import Dict

import rasterio as rio
from rasterio.errors import RasterioIOError
from rasterio.rpc import RPC

AVALAIBLE_FORMATS = ["geotiff", "geotiff_rpb", "rpb", "json"]


def rpc_writer(rpc_params: dict, filename: str, out_format: str = "geotiff", override: bool = False):
    """
    writer from  RPC dict, upper layer function to handle rpc_params for shareloc.rpc_reader.rpc_reader()
    :param rpc_params : rpc dict from rpc_reader
    :param filename : output file, in case of out_format is geotiff or geotiff_rpb, filename muse be in existing geotiff
    :param out_format : output_format in : geotiff, geotiff_rpb, rpb, json
    :param override : if True override rpc if already present in geotiff
    """
    if out_format not in AVALAIBLE_FORMATS:
        raise ValueError(f"{out_format} is no handled. must be one of {AVALAIBLE_FORMATS}")

    rio_rpcs = rpc_dict_to_rio_rpcs(rpc_params)

    if out_format.startswith("geotiff"):
        export_rpb = False
        if out_format.endswith("rpb"):
            export_rpb = True
        geotiff_rpc_updater(rio_rpcs, filename, export_rpb, override)
    elif out_format == "rpb":
        write_rio_rpc_as_rpb(rio_rpcs, filename)
    elif out_format == "json":
        write_rio_rpc_as_json(rio_rpcs, filename)


def geotiff_rpc_updater(rio_dict: Dict, filename, export_rpb=False, force_rpc=False):
    """
    geotiff updater

    :param rio_dict: path to geomodel
    :param filename: output filename
    :param export_rpb: write an .RPB file
    :param force_rpc: if False an exception is raised if filename already contains RPC
    """

    if export_rpb:
        rpb_path = os.path.splitext(filename)[0] + ".RPB"
        write_rio_rpc_as_rpb(rio_dict, rpb_path)
    else:
        try:
            with rio.open(filename, "r+") as src:
                rpcs = src.rpcs
                if rpcs is not None and force_rpc is False:
                    raise RuntimeError(f"{filename} already contain rpc")
                src.rpcs = RPC(**rio_dict)

        except RasterioIOError as rio_error:
            logging.debug("%s can not be read by rasterio", filename)
            logging.debug("    Rasterio error : %s", rio_error)


def write_rio_rpc_as_rpb(rpc: Dict, rpb_path: str, sat_id: str = "UNKNOWN", band_id: str = "P"):
    """
    Write .RPB file  (Rational Polynomial Coefficients)

    :param rpc: Geotiff rpc dict
    :param rpc_path: output file (.RPB)
    :param sat_id: sat ID (ex: "QB02")
    :param band_id: band ID (ex: "P")
    """

    def write_coeff(f, name: str, coeffs) -> None:
        """ "write coefficient bloc."""
        f.write(f"\t{name} = (\n")
        for i, c in enumerate(coeffs):
            sep = "," if i < len(coeffs) - 1 else ""
            f.write(f"\t\t{c}{sep}\n")
        f.write("\t);\n")

    with open(rpb_path, "w", encoding="utf-8") as f:
        f.write(f"satId =  {sat_id!r};\n")
        f.write(f"bandId = {band_id!r};\n")
        f.write('SpecId = "RPC00B";\n')
        f.write("BEGIN_GROUP = IMAGE\n")

        f.write("\terrBias = -1;\n")
        f.write("\terrRand = -1;\n")
        f.write(f"\tlineOffset = {rpc['line_off']};\n")
        f.write(f"\tsampOffset = {rpc['samp_off']};\n")
        f.write(f"\tlatOffset = {rpc['lat_off']};\n")
        f.write(f"\tlongOffset = {rpc['long_off']};\n")
        f.write(f"\theightOffset = {rpc['height_off']};\n")
        f.write(f"\tlineScale = {rpc['line_scale']};\n")
        f.write(f"\tsampScale = {rpc['samp_scale']};\n")
        f.write(f"\tlatScale = {rpc['lat_scale']};\n")
        f.write(f"\tlongScale = {rpc['long_scale']};\n")
        f.write(f"\theightScale = {rpc['height_scale']};\n")

        write_coeff(f, "lineNumCoef", rpc["line_num_coeff"])
        write_coeff(f, "lineDenCoef", rpc["line_den_coeff"])
        write_coeff(f, "sampNumCoef", rpc["samp_num_coeff"])
        write_coeff(f, "sampDenCoef", rpc["samp_den_coeff"])

        f.write("END_GROUP = IMAGE\n")
        f.write("END;\n")


def rpc_dict_to_rio_rpcs(rpc_params, topleftconvention=True) -> Dict:
    """
    Load via rasterio RPC object

    :param rpc_params: shareloc rpc reader dict
    :param topleftconvention: [0,0] position
    :type topleftconvention: boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)

    :return rio rpc dict
    """

    rio_rpc = {
        "height_off": rpc_params["offset_alt"],
        "height_scale": rpc_params["scale_alt"],
        "lat_off": rpc_params["offset_y"],
        "lat_scale": rpc_params["scale_y"],
        "line_den_coeff": rpc_params["den_row"],
        "line_num_coeff": rpc_params["num_row"],
        "line_off": rpc_params["offset_row"],
        "line_scale": rpc_params["scale_row"],
        "long_off": rpc_params["offset_x"],
        "long_scale": rpc_params["scale_x"],
        "samp_den_coeff": rpc_params["den_col"],
        "samp_num_coeff": rpc_params["num_col"],
        "samp_off": rpc_params["offset_col"],
        "samp_scale": rpc_params["scale_col"],
    }

    if topleftconvention:
        rio_rpc["samp_off"] -= 0.5
        rio_rpc["line_off"] -= 0.5

    return rio_rpc


def write_rio_rpc_as_json(rio_rpc_dict, filename):
    """
    Write rio rpc as json dict

    :param rio_rpc_dict: rio rpc dict RPC.to_dict()
    :param filename: json file
    """

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(rio_rpc_dict, f, ensure_ascii=False, indent=4)
