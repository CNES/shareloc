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

@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_eucl(col,row,alt):
    data_folder = test_path()

    fichier_direct_euclide = os.path.join(data_folder,'rpc/rpc_dir.euc')
    fichier_inverse_euclide =  os.path.join(data_folder,'rpc/rpc_inv.euc')

    fctrat = FonctRatD.from_euclidium(fichier_inverse_euclide, fichier_direct_euclide)

    (lon,lat, alt) = fctrat.direct_loc_h(row, col, alt)
    print (lon,lat)

    (row_ar,col_ar,__) = fctrat.inverse_loc(lon,lat,alt)
    print (col_ar,row_ar)
    assert(col_ar == pytest.approx(col,1e-2))
    assert(row_ar == pytest.approx(row,1e-2))

@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_phrdimap(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD.from_dimap(fichier_dimap)

    (lon,lat,alt) = fctrat.direct_loc_h(row,col,alt)
    print (lon,lat)

    (row_ar,col_ar,__) = fctrat.inverse_loc(lon,lat,alt)
    print (col_ar,row_ar)
    assert(col_ar == pytest.approx(col,1e-2))
    assert(row_ar == pytest.approx(row,1e-2))



@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_direct_inverse_iterative_vs_direct(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD.from_dimap(fichier_dimap)

    #(col,lig,alt)=(100,1000,400)
    (x0,y0, __) = fctrat.direct_loc_h(row,col,alt)
    (x_inv,y_inv) = fctrat.direct_loc_inverse_iterative(row,col,alt)
    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    assert (x0==pytest.approx(x_inv, 5.0/111111000))
    assert (y0==pytest.approx(y_inv, 5.0/111111000))



@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_direct_inverse_iterative(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD.from_dimap(fichier_dimap)

    (lon,lat,__) = fctrat.direct_loc_h(row,col,alt)
    (row_inv, col_inv,__) = fctrat.inverse_loc(lon, lat, alt)
    (lon_iter,lat_iter) = fctrat.direct_loc_inverse_iterative(row_inv,col_inv,alt)
    assert (lon==lon_iter)
    assert (lat==lat_iter)
