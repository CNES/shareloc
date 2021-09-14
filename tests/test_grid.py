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
Module to test functions that use direct grid
"""


import os
import pytest
import numpy as np
from utils import test_path

from shareloc.grid import Grid


@pytest.mark.unit_tests
def test_grid_eucl_vs_libgeorpc_drivers():
    """
    test grid readers
    """
    id_scene = "P1AP--2013072139303958CP"
    data_folder = test_path(alti="ellipsoide", id_scene=id_scene)
    gld = os.path.join(data_folder, "{}_H1.hd".format(id_scene))
    gri = Grid(gld)
    print(gri.nbrow)
    print(gri.nbcol)
    res_eucl = gri.direct_loc_h(100.0, 200.0, 100.0)

    data_folder_libgeo = os.path.join(os.environ["TESTPATH"], "ellipsoide", "{}_Libgeo".format(id_scene))
    gld_libgeo = os.path.join(data_folder_libgeo, "{}_H1.hdf".format(id_scene))
    print(gld_libgeo)
    gri_libgeo = Grid(gld_libgeo, grid_format="hdf")
    res_libgeo = gri_libgeo.direct_loc_h(100.0, 200.0, 100.0)
    np.testing.assert_allclose(res_eucl, res_libgeo, rtol=0, atol=1e-9)
