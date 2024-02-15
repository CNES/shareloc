#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
Shareloc Setup.py kept for compatibility and setuptools_scm configuration.
Main part is in setup.cfg file.
"""

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# For instance, to change compilation options:
# $> CC=gcc-11 CXX=g++-11 -march=native" pip install -e .[dev]
extensions = [
    Pybind11Extension(
        "bindings_cpp",
        [
            "shareloc/bindings/bind.cpp",
            "shareloc/bindings/dtm_intersection.cpp",
            "shareloc/bindings/rpc.cpp",
            "shareloc/bindings/GeoModelTemplate.cpp",
        ],
        language="c++",
        # cxx_std=20, # Uncomment to expect C++20
        extra_compile_args=["-O3"],
    )
]

setup(use_scm_version=True, cmdclass={"build_ext": build_ext}, ext_modules=extensions)
