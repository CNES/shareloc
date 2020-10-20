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

import os
import pytest
from utils import test_path

from shareloc.rpc.rpc_phr_v2 import FonctRatD

@pytest.mark.parametrize("col,lig,alt", [(600,200,125)])
def test_rpc_eucl(col,lig,alt):
    data_folder = test_path()

    fichier_direct_euclide = os.path.join(data_folder,'rpc/rpc_dir.euc')
    fichier_inverse_euclide =  os.path.join(data_folder,'rpc/rpc_inv.euc')

    fctrat = FonctRatD(fichier_direct_euclide, fichier_inverse_euclide)

    (lon,lat) = fctrat.evalue_loc_d(col,lig,alt)
    print (lon,lat)

    (col_ar,lig_ar) = fctrat.evalue_loc_i(lon,lat,alt)
    print (col_ar,lig_ar)
    assert(col_ar == pytest.approx(col,1e-2))
    assert(lig_ar == pytest.approx(lig,1e-2))

@pytest.mark.parametrize("col,lig,alt", [(600,200,125)])
def test_rpc_phrdimap(col,lig,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD(fichier_dimap)

    (lon,lat) = fctrat.evalue_loc_d(col,lig,alt)
    print (lon,lat)

    (col_ar,lig_ar) = fctrat.evalue_loc_i(lon,lat,alt)
    print (col_ar,lig_ar)
    assert(col_ar == pytest.approx(col,1e-2))
    assert(lig_ar == pytest.approx(lig,1e-2))


@pytest.mark.parametrize("col,lig,alt", [(600,200,125)])
def test_rpc_direct_inverse_iterative(col,lig,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD(fichier_dimap)

    #(col,lig,alt)=(100,1000,400)
    (x0,y0) = fctrat.evalue_loc_d(col,lig,alt)
    (x_inv,y_inv) = fctrat.evalue_loc_d_par_inversion_rpc_i(col,lig,alt)
    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    assert (x0==pytest.approx(x_inv, 5.0/111111000))
    assert (y0==pytest.approx(y_inv, 5.0/111111000))
