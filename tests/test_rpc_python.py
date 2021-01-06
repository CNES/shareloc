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
import numpy as np

from shareloc.rpc.rpc import RPC, identify_dimap, identify_ossim_kwl



def test_rpc_drivers():
    data_folder = test_path()

    file_direct_euclide = os.path.join(data_folder,'rpc/rpc_dir.euc')
    file_inverse_euclide =  os.path.join(data_folder,'rpc/rpc_inv.euc')

    fctrat_eucl = RPC.from_any(file_inverse_euclide, file_direct_euclide)

    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat_dimap = RPC.from_any(file_dimap)

    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_geom = os.path.join(data_folder,'rpc/{}.geom'.format(id_scene))
    fctrat_geom = RPC.from_any(file_geom)

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

@pytest.mark.parametrize("id_scene,lon,lat,alt, col_vt,row_vt", [('PHR1B_P_201709281038393_SEN_PRG_FC_178609-001',
                                                                      7.048662660737769592, 43.72774839443545858, 0.0,
                                                                      100.0, 200.0),
                                        ('PHR1B_P_201709281038045_SEN_PRG_FC_178608-001',
                                         7.11526088296757386331137240632,43.684281179313565246502548689,850.0,10121.0657,10235.9310)])
def test_rpc_ossim_kwl(id_scene,lon,lat,alt,row_vt,col_vt):
    data_folder = test_path()
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    fctrat_geom = RPC.from_any(file_geom, topleftconvention=True)
    print("{} {} {}".format(lon,lat,alt))
    (row,col,__) = fctrat_geom.inverse_loc(lon,lat,alt)
    print("col {} row {}".format(col,row))
    assert(col == pytest.approx(col_vt, abs = 1e-2))
    assert(row == pytest.approx(row_vt, abs = 1e-2))



@pytest.mark.parametrize("lon,lat,alt", [(7.048662660737769592,43.72774839443545858,0.0)])
def test_rpc_eucl(lon,lat,alt):
    data_folder = test_path()
    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_inverse = os.path.join(data_folder, 'rpc/{}_inv.txt'.format(id_scene))
    fctrat_eucl = RPC.from_euclidium(file_inverse)
    print("{} {} {}".format(lon,lat,alt))
    (row,col,__) = fctrat_eucl.inverse_loc(lon,lat,alt)
    print("col {} row {}".format(col,row))
    assert(col == pytest.approx(100.5, abs = 1e-2))
    assert(row == pytest.approx(200.5, abs = 1e-2))


@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_phrdimap(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

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

    fctrat = RPC.from_dimap_v1(file_dimap)

    #(col,lig,alt)=(100,1000,400)
    (x0,y0, __) = fctrat.direct_loc_h(row,col,alt)
    (x_inv,y_inv,__) = fctrat.direct_loc_inverse_iterative(row,col,alt)
    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    assert (x0==pytest.approx(x_inv, abs = 10.0/111111000))
    assert (y0==pytest.approx(y_inv, abs =10.0/111111000))


def test_rpc_direct_inverse_iterative_vs_direct_multiple_points():
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (col,row,alt)=(np.array([600, 610]), np.array([200, 210]), np.array([125]))
    P0 = fctrat.direct_loc_h(row,col,alt)
    (rowinv,colinv,__) = fctrat.inverse_loc(P0[:,0],P0[:,1],alt)
    P1 = fctrat.direct_loc_inverse_iterative(rowinv,colinv,alt)

    #print("comparaison loc directe RPC / loc inverse RPC iterative""")
    #"""Les erreurs constates sont dus aux differences entre loc dir et loc inv RPC"""
    #Point error col = 600, row = 200
    assert (P0[0, 0]==P1[0][0])
    assert (P0[0, 1]==P1[1][0])

    # Point error col = 601, row = 201
    assert (P0[1, 0]==P1[0][1])
    assert (P0[1, 1]==P1[1][1])


def test_rpc_direct_iterative_nan():
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (col,row,alt)=(np.array([600, np.nan]), np.array([200, 210]), np.array([125]))
    P0 = fctrat.direct_loc_inverse_iterative(row,col,alt, fill_nan = True)
    P1 = fctrat.direct_loc_inverse_iterative(row[0], col[0], alt,fill_nan = True)
    #Point error col = 600, row = 200
    assert (P0[0][0]==P1[0])
    assert (P0[0][1]==P1[1])

def test_rpc_direct_iterative_nan():
    # comparasion between tabs en tab containing Nan cf. issue #46
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)


    (col,row,alt)=(np.array([600, np.nan]), np.array([200, 210]), np.array([125]))
    direct_loc_tab = fctrat.direct_loc_inverse_iterative(row,col,alt, fill_nan = True)
    direct_loc_one_value = fctrat.direct_loc_inverse_iterative(row[0], col[0], alt, fill_nan = True)
    (col,row,alt)=(np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), np.array([125]))
    direct_loc_nan = fctrat.direct_loc_inverse_iterative(row,col,alt,10,False)
    assert(np.all(np.isnan(direct_loc_nan[0])))
    assert(np.all(np.isnan(direct_loc_nan[1])))

    #Point error col = 600, row = 200
    assert (direct_loc_tab[0][0]==direct_loc_one_value[0][0])
    assert (direct_loc_tab[1][0]==direct_loc_one_value[1][0])
    assert (direct_loc_tab[0][1]==fctrat.offset_X)
    assert (direct_loc_tab[1][1]==fctrat.offset_Y)


@pytest.mark.parametrize("col,row,alt", [(600,200,125)])
def test_rpc_direct_inverse_iterative(col,row,alt):
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    file_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = RPC.from_dimap_v1(file_dimap)

    (lon,lat,__) = fctrat.direct_loc_h(row,col,alt)
    (row_inv, col_inv,__) = fctrat.inverse_loc(lon, lat, alt)
    (lon_iter,lat_iter,__) = fctrat.direct_loc_inverse_iterative(row_inv,col_inv,alt)
    assert (lon==lon_iter)
    assert (lat==lat_iter)

def test_rpc_minmax():
    data_folder = test_path()
    id_scene = 'P1BP--2018122638935449CP'
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))
    fctrat = RPC.from_any(fichier_dimap)
    (h_min,h_max) = fctrat.get_alt_min_max()
    assert (h_min == 532.5)
    assert (h_max == 617.5)
