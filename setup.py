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
import os

# from setuptools.command.build_ext import build_ext
from pybind11.setup_helpers import Pybind11Extension, build_ext, intree_extensions
from setuptools import Extension, setup

# create libs folder to contain .so files (if it doesn't exist)
if not os.path.exists("libs"):
    os.makedirs("libs")

# Main setup with setup.cfg file.
extensions = [Pybind11Extension("libs.pbrpc", ["shareloc/draft_pybind/bind_helloworld.cpp"])]  # "lib.pbhelloworld"

# extensions = intree_extensions("pbrpc", ["shareloc/draft_pybind/hello_world.cpp"])

setup(use_scm_version=True, cmdclass={"build_ext": build_ext}, ext_modules=extensions)
