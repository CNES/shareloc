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
shareloc-rpc-converter
"""

import argparse
import logging
import textwrap

from shareloc.geomodels.rpc_readers import rpc_reader
from shareloc.geomodels.rpc_writers import rpc_writer


def rpc_converter(input_rpc: str, output_file: str, rpc_format: str, loglevel):
    """
    convert rpc to another format.
    :param input_rpc: input rpc
    :param output_file: output file, existing geotiff or filename
    :param rpc_format: output rpc format
    """
    setup_logging(loglevel)
    rpc_dict = rpc_reader(input_rpc)
    driver_type = rpc_dict["driver_type"]
    logging.info("driver type used to read rpc %s", driver_type)

    rpc_writer(rpc_dict, output_file, rpc_format)


def setup_logging(loglevel=logging.INFO):
    logging.getLogger().setLevel(loglevel)


def cli():
    """
    Command Line Interface
    """

    parser = argparse.ArgumentParser(
        "shareloc-rpcconverter",
        description="convert RPC models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
    This script takes as input an RPC model and transforms it to another format"
    """
        ),
    )
    parser.add_argument("input_rpc", type=str, help="input rpc")
    parser.add_argument("output_file", type=str, help="output_rfile")
    parser.add_argument(
        "rpc_format",
        type=str,
        help="output rpc format in geotiff, geotiff_rpb, rpb, json. "
        "if format is geotiff or geotiff_rpb output_file "
        "must be an existing tif, which will be updated with rpc tifftag,"
        "or external .RPB file in case of geotiff_rpb, if rpb, output_file"
        " must be a .RPB.",
    )
    parser.add_argument(
        "--loglevel",
        default="DEBUG",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Logger level (default: WARNING. Should be one of " "(DEBUG, INFO, WARNING, ERROR, CRITICAL)",
    )

    args = parser.parse_args()

    rpc_converter(**vars(args))


if __name__ == "__main__":
    cli()
