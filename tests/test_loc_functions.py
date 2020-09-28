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
import numpy as np
from utils import test_path

from shareloc.grid import grid, coloc
from shareloc.dtm import DTM
from shareloc.rpc.rpc_phr_v2 import FonctRatD
from shareloc.localization import Localization


def prepare_loc(alti = 'geoide', id_scene='P1BP--2017030824934340CP'):
    """
    Read multiH grid

    :param alti: alti validation dir
    :param id_scene: scene euclidium ID
    :return: multi H grid
    :rtype: str
    """   
    data_folder = test_path(alti, id_scene)
    
    #chargement du mnt
    fic = os.path.join(data_folder,'MNT_extrait/mnt_extrait.c1')
    dtmbsq = DTM(fic)
    
    #chargement des grilles
    gld = os.path.join(data_folder,'grilles_gld_xH/{}_H1.hd'.format(id_scene))
    gri = grid(gld)
    
   
    return dtmbsq,gri

"""
calcul_gld3d -all -path_visee ./grilles_gld_xH -nom_visee P1BP--2017030824934340CP_H1 -type_visee LocalisateurGrille_Directe \
-mnt ./MNT_extrait/mnt_extrait -repter GRS80:G-D/:H-M -path_grille . -nom_grille test_intersect_euclide -convention BABEL \
-format BSQ -nbcol 200 -nblig 200 -pascol 50 -paslig 60 -j0 0 -i0 0 -col0 100.5 -lig0 100.5 -matrice 1 0 0 1
"""


@pytest.mark.parametrize("idxcol,idxlig", [(20,10)])
@pytest.mark.parametrize("lig0,col0,paslig,pascol,nblig,nbcol", [(100.5,150.5,60,50,200,200)])
@pytest.mark.unit_tests
def test_gld_dtm(idxlig,idxcol,lig0,col0,paslig,pascol,nblig,nbcol):
    """
    Test loc direct grid on dtm function
    """
    dtmbsq,gri = prepare_loc()

    gri_gld = gri.direct_loc_grid_dtm(lig0, col0, paslig, pascol, nblig, nbcol, dtmbsq)

    lonlatalt = gri_gld[:,idxlig,idxcol]

    lig = lig0 + paslig * idxlig
    col = col0 + pascol * idxcol

    valid_lonlatalt = gri.direct_loc_dtm(lig, col, dtmbsq)
    print("lon {} lat {} alt {} ".format(lonlatalt[0], lonlatalt[1], lonlatalt[2]))
    print("valid lon {} lat {} alt {} ".format(valid_lonlatalt[0], valid_lonlatalt[1], valid_lonlatalt[2]))
    assert(lonlatalt == pytest.approx(valid_lonlatalt, abs=1e-12))


@pytest.mark.parametrize("col,lig", [(50.5,100.5)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.21700176041541,21.959197148974,238.0)])
@pytest.mark.unit_tests
def test_loc_dir_check_cube_dtm(col,lig,valid_lon,valid_lat,valid_alt):
    """
    Test direct localization check dtm cube
    """
    dtmbsq,gri = prepare_loc()

    visee = np.zeros((3, gri.nbalt))
    vislonlat = gri.interpolate_grid_in_plani(lig, col)
    visee[0,:] = vislonlat[0]
    visee[1,:] = vislonlat[1]
    visee[2,:] = gri.alts_down
    v = visee.T
    (code1, code2, PointB, dH3D) = dtmbsq.checkCubeDTM(v)
    assert(PointB[0] == pytest.approx(valid_lon,abs=1e-12))
    assert(PointB[1] == pytest.approx(valid_lat,abs=1e-12))
    assert(PointB[2] == pytest.approx(valid_alt,abs=1e-12))

@pytest.mark.parametrize("col,lig", [(50.5,100.5)])
@pytest.mark.parametrize("valid_lon,valid_lat", [(57.2169100597702,21.96277930762832)])
@pytest.mark.unit_tests
def test_loc_dir_interp_visee_unitaire_gld(lig,col,valid_lon,valid_lat):
    """
    Test los interpolation
    """
    ___, gri = prepare_loc()
    visee = gri.interpolate_grid_in_plani(lig, col)
    print(visee)
    assert (visee[0][1] == pytest.approx(valid_lon,abs=1e-12))
    assert (visee[1][1] == pytest.approx(valid_lat, abs=1e-12))



@pytest.mark.parametrize("col,lig,h", [(50.5,100.5,100.0)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.2170054518422,21.9590529453258,100.0)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_h(col,lig,h,valid_lon,valid_lat,valid_alt):
    """
    Test direct localization at constant altitude
    """
    ___,gri = prepare_loc()
    loc = Localization(gri, dtm = None)
    lonlatalt = loc.direct(lig, col, h)

    diff_lon = lonlatalt[0] - valid_lon
    diff_lat = lonlatalt[1] - valid_lat
    diff_alt = lonlatalt[2] - valid_alt
    print("direct localization at constant altitude lig : {} col {} alt {}".format(lig,col,h))
    print("lon {} lat {} alt {} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))
    print('diff_lon {} diff_lat {} diff_alt {}'.format(diff_lon, diff_lat, diff_alt))
    assert(valid_lon == pytest.approx(lonlatalt[0],abs=1e-12))
    assert(valid_lat == pytest.approx(lonlatalt[1],abs=1e-12))
    assert(valid_alt == pytest.approx(lonlatalt[2],abs=1e-8))


@pytest.mark.parametrize("col,lig,h", [(150.5,100,100.0)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_vs_loc_rpc(lig, col, h):
    """
    Test direct localization coherence with grid and RPC
    """
    id_scene = 'P1BP--2018122638935449CP'
    ___,gri = prepare_loc('ellipsoide',id_scene)
    loc_grid = Localization(gri)
    #init des predicteurs
    lonlatalt = loc_grid.direct(lig, col, h)

    data_folder = test_path()
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD(fichier_dimap)

    loc_rpc = Localization(fctrat)
    lonlatalt_rpc = loc_rpc.direct(lig + 0.5, col + 0.5, h)

    diff_lon = lonlatalt[0] - lonlatalt_rpc[0]
    diff_lat = lonlatalt[1] - lonlatalt_rpc[1]
    diff_alt = lonlatalt[2] - lonlatalt_rpc[2]
    assert(diff_lon == pytest.approx(0.0, abs=1e-7))
    assert(diff_lat == pytest.approx(0.0, abs=1e-7))
    assert(diff_alt == pytest.approx(0.0, abs=1e-7))


@pytest.mark.parametrize("index_x,index_y", [(10.5,20.5)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_dtm(index_x,index_y):
    """
    Test direct localization on DTM
    """
    dtmbsq, gri = prepare_loc()
    loc = Localization(gri, dtm=dtmbsq)


    vect_index = [index_x, index_y]
    [lon,lat] = dtmbsq.DTMToTer(vect_index)
    print([lon,lat])
    alt = dtmbsq.Interpolate(index_x - 0.5, index_y - 0.5)

    lig, col, valid = loc.inverse(lon, lat, alt)
    print("lig col ",lig,col)
    lonlath = loc.direct(lig, col)
    assert(lon == pytest.approx(lonlath[0],abs=1e-8))
    assert(lat == pytest.approx(lonlath[1],abs=1e-8))
    assert(alt == pytest.approx(lonlath[2],abs=1e-4))

@pytest.mark.parametrize("col,lig,h", [(50.5,100.5,100.0)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.2170054518422,21.9590529453258,100.0)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_h(col,lig,h,valid_lon,valid_lat,valid_alt):
    """
    Test direct localization at constant altitude
    """
    ___,gri = prepare_loc()
    loc = Localization(gri, dtm = None)
    lonlatalt = loc.direct(lig, col, h)

    diff_lon = lonlatalt[0] - valid_lon
    diff_lat = lonlatalt[1] - valid_lat
    diff_alt = lonlatalt[2] - valid_alt
    print("direct localization at constant altitude lig : {} col {} alt {}".format(lig,col,h))
    print("lon {} lat {} alt {} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))
    print('diff_lon {} diff_lat {} diff_alt {}'.format(diff_lon, diff_lat, diff_alt))
    assert(valid_lon == pytest.approx(lonlatalt[0],abs=1e-12))
    assert(valid_lat == pytest.approx(lonlatalt[1],abs=1e-12))
    assert(valid_alt == pytest.approx(lonlatalt[2],abs=1e-8))


@pytest.mark.parametrize("valid_lig,valid_col", [(50.5,10.0)])
@pytest.mark.parametrize("lon,lat,alt", [(57.2167252772905,21.9587514585812,10.0)])
@pytest.mark.unit_tests
def test_sensor_loc_inv(lon,lat,alt,valid_col,valid_lig):
    """
    Test inverse localization
    """

    ___,gri = prepare_loc()

    loc = Localization(gri)
    inv_lig,inv_col,valid = loc.inverse(lon,lat,alt)

    print("inverse localization  : lon {} lat {} alt {}".format(lon,lat,alt))
    print("lig {} col {}  ".format(inv_lig, inv_col))
    print('diff_lig {} diff_col {} '.format(inv_lig - valid_lig, inv_col - valid_col))   
    assert(inv_lig == pytest.approx(valid_lig,abs=1e-2))
    assert(inv_col == pytest.approx(valid_col,abs=1e-2))


@pytest.mark.parametrize("lon,lat,alt", [(2.12026631, 31.11245154,10.0)])
@pytest.mark.unit_tests
def test_sensor_loc_inv_vs_loc_rpc(lon, lat, alt):
    """
    Test direct localization coherence with grid and RPC
    """
    id_scene = 'P1BP--2018122638935449CP'
    ___,gri = prepare_loc('ellipsoide',id_scene)
    loc_grid = Localization(gri)
    #init des predicteurs
    [row, col, valid] = loc_grid.inverse(lon, lat, alt)
    data_folder = test_path()
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD(fichier_dimap)

    loc_rpc = Localization(fctrat)
    [row_rpc, col_rpc, valid] = loc_rpc.inverse(lon, lat, alt)
    diff_row = row_rpc - row
    diff_col = col_rpc - col
    assert (diff_row == pytest.approx(0.5, abs=1e-2))
    assert (diff_col == pytest.approx(0.5, abs=1e-2))




@pytest.mark.parametrize("col_min_valid", [([4.15161251e-02, 1.95057636e-01, 1.10977819e+00, -8.35016563e-04,-3.50772271e-02, -9.46432481e-03])])
@pytest.mark.parametrize("lig_min_valid", [([ 0.05440845,  1.26513831, -0.36737151, -0.00229532, -0.07459378, -0.02558954])])
@pytest.mark.parametrize("col_max_valid", [([1.76451389e-02, 2.05533045e-01, 1.11758291e+00, -9.50086076e-04,-3.59923603e-02, -1.03291594e-02])])
@pytest.mark.parametrize("lig_max_valid", [([0.07565692, 1.27499912, -0.36677813, -0.00252395, -0.07539624, -0.0270914])])
@pytest.mark.parametrize("valid_offset_lon",[([57.37295223744326, 0.15660032225072484])])
@pytest.mark.parametrize("valid_offset_lat",[([22.066877016445275, 0.14641205050748773])])
@pytest.mark.parametrize("valid_offset_lig",[([24913.0, 24912.5])])
@pytest.mark.parametrize("valid_offset_col",[([19975.5, 19975.0])])
@pytest.mark.unit_tests
def test_pred_loc_inv(col_min_valid, lig_min_valid,col_max_valid,lig_max_valid,valid_offset_lon,valid_offset_lat,valid_offset_lig,valid_offset_col):
    """
    Test inverse localization
    """
    #init des predicteurs
    ___,gri = prepare_loc()
    gri.estimate_inverse_loc_predictor()

    assert (gri.pred_col_min.flatten() == pytest.approx(col_min_valid, abs=1e-6))
    assert (gri.pred_lig_min.flatten() == pytest.approx(lig_min_valid, abs=1e-6))
    assert (gri.pred_col_max.flatten() == pytest.approx(col_max_valid, abs=1e-6))
    assert (gri.pred_lig_max.flatten() == pytest.approx(lig_max_valid, abs=1e-6))
    assert(gri.pred_ofset_scale_lon == pytest.approx(valid_offset_lon, abs=1e-12))
    assert(gri.pred_ofset_scale_lat == pytest.approx(valid_offset_lat, abs=1e-12))
    assert(gri.pred_ofset_scale_lig == pytest.approx(valid_offset_lig, abs=1e-6))
    assert(gri.pred_ofset_scale_col == pytest.approx(valid_offset_col, abs=1e-6))

@pytest.mark.parametrize("col,lig", [(50.5,100.5)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.21700367698209,21.95912227930429,166.351227229112)])
@pytest.mark.unit_tests
def test_loc_intersection(lig,col,valid_lon,valid_lat,valid_alt):
    """
    Test direct localization intersection function
    """
    dtmbsq,gri = prepare_loc()

    visee = np.zeros((3, gri.nbalt))
    vislonlat = gri.interpolate_grid_in_plani(lig, col)
    visee[0,:] = vislonlat[0]
    visee[1,:] = vislonlat[1]
    visee[2,:] = gri.alts_down
    v = visee.T
    (code1, code2, PointB, dH3D) = dtmbsq.checkCubeDTM(v)
    (code3,code4,Point_dtm) = dtmbsq.intersection(v, PointB, dH3D)
    assert(Point_dtm[0] == pytest.approx(valid_lon,abs=1e-12))
    assert(Point_dtm[1] == pytest.approx(valid_lat,abs=1e-12))
    assert(Point_dtm[2] == pytest.approx(valid_alt,abs=1e-12))



@pytest.mark.parametrize("col,lig,h", [(20.5,150.5,10.0)])
@pytest.mark.unit_tests
def test_loc_dir_loc_inv(lig, col, h):
    """
    Test direct localization followed by inverse one
    """
    ___,gri = prepare_loc()
    #init des predicteurs
    gri.estimate_inverse_loc_predictor()
    (lon,lat,alt) = gri.direct_loc_h(lig, col, h)
    inv_lig,inv_col,valid = gri.inverse_loc(lon,lat,alt)

    print('lig {} col {} valid {}'.format(inv_lig, inv_col, valid))
    assert(lig == pytest.approx(inv_lig,abs=1e-2))
    assert(col == pytest.approx(inv_col,abs=1e-2))
    assert(valid == 1)

@pytest.mark.parametrize("col,lig,h", [(150.5,20.5,10.0)])
@pytest.mark.unit_tests
def test_loc_dir_loc_inv_rpc(lig, col, h):
    """
    Test direct localization followed by inverse one
    """
    id_scene = 'P1BP--2018122638935449CP'
    ___,gri = prepare_loc('ellipsoide',id_scene)
    #init des predicteurs
    gri.estimate_inverse_loc_predictor()
    lonlatalt = gri.direct_loc_h(lig, col, h)

    data_folder = test_path()
    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP_{}.XML'.format(id_scene))

    fctrat = FonctRatD(fichier_dimap)
    (inv_lig, inv_col,__) = fctrat.inverse_loc(lonlatalt[0], lonlatalt[1], lonlatalt[2])
    print('lig {} col {}'.format(inv_lig, inv_col))

    assert(lig == pytest.approx(inv_lig - 0.5, abs=1e-2))
    assert(col == pytest.approx(inv_col - 0.5, abs=1e-2))


@pytest.mark.parametrize("l0_src,c0_src, paslig_src, pascol_src,nblig_src,nbcol_src", [(0.5,1.5,10,100,20, 20)])
@pytest.mark.parametrize("col,lig", [(1,3)])
@pytest.mark.unit_tests
def test_coloc(l0_src, c0_src,paslig_src,pascol_src,nblig_src,nbcol_src,lig,col):
    """
    Test coloc function
    """
    dtmbsq,gri = prepare_loc()
    gri.estimate_inverse_loc_predictor()

    gricol = coloc(gri, gri, dtmbsq, l0_src, c0_src, paslig_src, pascol_src, nblig_src, nbcol_src)

    assert(gricol[0, lig, col] == pytest.approx(lig * paslig_src + l0_src,1e-6))
    assert(gricol[1, lig, col] == pytest.approx(col * pascol_src + c0_src, 1e-6))

@pytest.mark.parametrize("index_x,index_y", [(10.5,20.5)])
@pytest.mark.parametrize("valid_alt", [(198.0)])
@pytest.mark.unit_tests
def test_interp_dtm(index_x,index_y,valid_alt):
    """
    Test coloc function
    """
    dtmbsq, ___ = prepare_loc()

    vect_index = [index_x, index_y]
    coords = dtmbsq.DTMToTer(vect_index)
    print(coords)
    alti = dtmbsq.Interpolate(index_x - 0.5, index_y - 0.5)
    assert(alti == valid_alt)


    
    



