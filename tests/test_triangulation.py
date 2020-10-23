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
import xarray as xr

import time

from shareloc.grid import grid, fct_coloc
from shareloc.dtm import DTM
from shareloc.triangulation.triangulation import sensor_triangulation, epipolar_triangulation, transform_disp_to_matches
from shareloc.rectification.rectification_grid import rectification_grid

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

    
@pytest.mark.parametrize("col,lig,h", [(1000.5,1500.5,10.0)])
@pytest.mark.unit_tests
def test_sensor_triangulation(lig, col, h):
    """
    Test sensor triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___,gri_right = prepare_loc('ellipsoide',id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc('ellipsoide',id_scene_left)
    #init des predicteurs
    gri_right.init_pred_loc_inv()
    lonlatalt = gri_left.fct_locdir_h(lig, col, h)
    inv_lig,inv_col,valid = gri_right.fct_locinv(lonlatalt)

    print('lig {} col {} valid {}'.format(inv_lig, inv_col, valid))
    matches = np.zeros([1,4])
    matches[0,:] = [lig,col,inv_lig,inv_col]
    #matches[1,:] = [lig + 10, col + 5, inv_lig + 12, inv_col + 7]

    point_ecef, point_wgs84 = sensor_triangulation(matches,gri_left,gri_right)
    assert(lonlatalt[0] == pytest.approx(point_wgs84[0,0],abs=1e-8))
    assert(lonlatalt[1] == pytest.approx(point_wgs84[0,1],abs=1e-8))
    assert(lonlatalt[2] == pytest.approx(point_wgs84[0,2],abs=6e-3))
    #assert(col == pytest.approx(inv_col,abs=1e-2))
    #assert(valid == 1)




@pytest.mark.unit_tests
def test_epi_triangulation_sift():
    """
    Test epipolar triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___,gri_right = prepare_loc('ellipsoide',id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc('ellipsoide',id_scene_left)



    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_left))
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))


    matches_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches.npy")
    matches = np.load(matches_filename)

    print("matches {}".format(matches[1:10,0:2]))

    point_ecef, point_wgs84 = epipolar_triangulation(matches, None,'sift',gri_left,gri_right,grid_left_filename,grid_right_filename)

    print(point_wgs84[1:10,:])

    #assert(lonlatalt[0] == pytest.approx(point_wgs84[0,0],abs=1e-8))
    #assert(lonlatalt[1] == pytest.approx(point_wgs84[0,1],abs=1e-8))
    ##assert(lonlatalt[2] == pytest.approx(point_wgs84[0,2],abs=6e-3))
    #assert(col == pytest.approx(inv_col,abs=1e-2))
    #assert(valid == 1)



@pytest.mark.unit_tests
def test_epi_triangulation_disp():
    """
     Test epipolar triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___, gri_right = prepare_loc('ellipsoide', id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc('ellipsoide', id_scene_left)

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
                                      "grid_{}.tif".format(id_scene_left))
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
                                       "grid_{}.tif".format(id_scene_right))


    disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)

    point_ecef, point_wgs84 = epipolar_triangulation(disp, None, 'disp', gri_left, gri_right, grid_left_filename,
                                                   grid_right_filename)
    array_shape = disp.disp.values.shape
    array_epi_wgs84 = point_wgs84.reshape((array_shape[0], array_shape[1],3))
    pc_dataset = xr.Dataset({'pc_wgs84': (['row', 'col','coords'], array_epi_wgs84)},
                          coords={'row': disp.coords['row'], 'col': disp.coords['col'], 'coords' : np.arange(3).astype(np.uint8)})
    disp = xr.merge((disp, pc_dataset))
    #point_ecef, point_wgs84 = epipolar_triangulation(matches, None,'sift',gri_left,gri_right,grid_left_filename,grid_right_filename)

    #print(point_wgs84[1:10,:])

    #assert(lonlatalt[0] == pytest.approx(point_wgs84[0,0],abs=1e-8))
    #assert(lonlatalt[1] == pytest.approx(point_wgs84[0,1],abs=1e-8))
    ##assert(lonlatalt[2] == pytest.approx(point_wgs84[0,2],abs=6e-3))
    #assert(col == pytest.approx(inv_col,abs=1e-2))
    #assert(valid == 1)


