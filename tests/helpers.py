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
module for test utilities
"""

import os


def data_path(alti="", scene_id=""):
    """
    return the data path, when used without any argument data_path() returns data directory
    :param alti: first sub dir corresponding to datum ("ellipsoide" or "geoid")
    :type alti: str
    :param scene_id: second sub dir corresponding to the scene id
    :type scene_id: str
    :return: data path.
    :rtype: str
    """
    data_root_folder = os.path.join(os.path.dirname(__file__), "data")
    sub_folder = os.path.join(alti, scene_id)
    return os.path.join(data_root_folder, sub_folder)
