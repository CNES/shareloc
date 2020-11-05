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

from shareloc.rpc.rpc_phr_v2 import FonctRatD, identify_dimap, identify_ossim_kwl



def test_rpc_drivers():
    data_folder = test_path()

    file_direct_euclide = os.path.join(data_folder,'rpc/rpc_dir.euc')
    file_inverse_euclide =  os.path.join(data_folder,'rpc/rpc_inv.euc')

    fctrat_eucl = FonctRatD.from_any(file_inverse_euclide, file_direct_euclide)

    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat_dimap = FonctRatD.from_any(file_dimap)

    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_geom = os.path.join(data_folder,'rpc/{}.geom'.format(id_scene))
    fctrat_geom = FonctRatD.from_any(file_geom)

    assert (fctrat_eucl.driver_type == 'euclidium')
    assert (fctrat_dimap.driver_type == 'dimap_v1.4')
    assert (fctrat_geom.driver_type == 'ossim_kwl')

def test_identify_ossim_kwl():
    data_folder = test_path()
    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    ossim_model =  identify_ossim_kwl(file_geom)
    assert (ossim_model == 'ossimPleiadesModel')



def test_identify_dimap():
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder, 'rpc/PHRDIMAP_{}.XML'.format(id_scene))
    dimap_version =  identify_dimap(file_dimap)
    assert (dimap_version == '1.4')

@pytest.mark.parametrize("lon,lat,alt", [(7.048662660737769592,43.72774839443545858,0.0)])
def test_rpc_ossim_kwl(lon,lat,alt):
    data_folder = test_path()
    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    fctrat_geom = FonctRatD.from_any(file_geom)
    print("{} {} {}".format(lon,lat,alt))
    (row,col,__) = fctrat_geom.inverse_loc(lon,lat,alt)
    print("col {} row {}".format(col,row))
    assert(col == pytest.approx(100, abs = 1e-2))
    assert(row == pytest.approx(200, abs = 1e-2))



@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_eucl(col,row,alt):
    data_folder = test_path()

    file_direct_euclide = os.path.join(data_folder,'rpc/rpc_dir.euc')
    file_inverse_euclide =  os.path.join(data_folder,'rpc/rpc_inv.euc')

    fctrat = FonctRatD.from_euclidium(file_inverse_euclide, file_direct_euclide)

    (lon,lat, alt) = fctrat.direct_loc_h(row, col, alt)
    print (lon,lat)

    (row_ar,col_ar,__) = fctrat.inverse_loc(lon,lat,alt)
    print (col_ar,row_ar)
    assert(col_ar == pytest.approx(col,abs = 1e-2))
    assert(row_ar == pytest.approx(row,abs = 1e-2))

@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_phrdimap(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD.from_dimap_v1(file_dimap)

    (lon,lat,alt) = fctrat.direct_loc_h(row,col,alt)
    print (lon,lat)

    (row_ar,col_ar,__) = fctrat.inverse_loc(lon,lat,alt)
    print (col_ar,row_ar)
    assert(col_ar == pytest.approx(col,abs = 2e-2))
    assert(row_ar == pytest.approx(row,abs = 2e-2))



@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_direct_inverse_iterative_vs_direct(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD.from_dimap_v1(file_dimap)

    #(col,lig,alt)=(100,1000,400)
    (x0,y0, __) = fctrat.direct_loc_h(row,col,alt)
    (x_inv,y_inv) = fctrat.direct_loc_inverse_iterative(row,col,alt)
    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    assert (x0==pytest.approx(x_inv, abs = 10.0/111111000))
    assert (y0==pytest.approx(y_inv, abs =10.0/111111000))



@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_direct_inverse_iterative(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD.from_dimap_v1(file_dimap)

    (lon,lat,__) = fctrat.direct_loc_h(row,col,alt)
    (row_inv, col_inv,__) = fctrat.inverse_loc(lon, lat, alt)
    (lon_iter,lat_iter) = fctrat.direct_loc_inverse_iterative(row_inv,col_inv,alt)
    assert (lon==lon_iter)
    assert (lat==lat_iter)
