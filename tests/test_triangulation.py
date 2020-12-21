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
from shareloc.localization import Localization
from shareloc.grid import grid
from shareloc.dtm import DTM
from shareloc.triangulation.triangulation import distance_point_los,sensor_triangulation
from shareloc.triangulation.triangulation import epipolar_triangulation, transform_disp_to_matches
from shareloc.rectification.rectification_grid import rectification_grid
from shareloc.rpc.rpc import RPC

import rasterio as rio


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

    
@pytest.mark.parametrize("col,row,h", [(1000.5,1500.5,10.0)])
@pytest.mark.unit_tests
def test_sensor_triangulation(row, col, h):
    """
    Test sensor triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___,gri_right = prepare_loc('ellipsoide',id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc('ellipsoide',id_scene_left)
    #init des predicteurs
    gri_right.estimate_inverse_loc_predictor()
    lonlatalt = gri_left.direct_loc_h(row, col, h)

    inv_row,inv_col,valid = gri_right.inverse_loc(lonlatalt[0],lonlatalt[1],lonlatalt[2])


    matches = np.zeros([1,4])
    matches[0,:] = [col,row,inv_col,inv_row]
    #matches[1,:] = [lig + 10, col + 5, inv_lig + 12, inv_col + 7]

    point_ecef, point_wgs84, distance  = sensor_triangulation(matches,gri_left,gri_right, residues = True)

    assert(lonlatalt[0] == pytest.approx(point_wgs84[0,0],abs=1e-8))
    assert(lonlatalt[1] == pytest.approx(point_wgs84[0,1],abs=1e-8))
    assert(lonlatalt[2] == pytest.approx(point_wgs84[0,2],abs=6e-3))
    assert(distance == pytest.approx(0.0,abs=1e-3))
    #assert(valid == 1)

@pytest.mark.unit_tests
def test_triangulation_residues():

    class simulatedLOS:
        """ line of sight class
        """

        def __init__(self):
            self.sis=np.array([[100.0, 10.0, 200.0],[100.0, 10.0, 200.0]])
            self.vis=np.array([[0.0, 1., 0.0],[0.0, 1., 0.0] ])
    los = simulatedLOS()

    distance = 10.0
    point = los.sis + 100.0 * los.vis + distance* np.array([[0.0, 0.0, 1.0],[0.0, 0.0, 1.0]])
    residue = distance_point_los(los, point)
    assert (distance == pytest.approx(residue, abs=1e-9))

@pytest.mark.unit_tests
def test_epi_triangulation_sift():
    """
    Test epipolar triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___,gri_right = prepare_loc('ellipsoide',id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc('ellipsoide',id_scene_left)

    #row = 10.0
    #col = 1000.0
    #h = 10.0
    #loc_left = Localization(grid=gri_left, dtm=None)
    #lonlatalt = loc_left.forward(row, col, h)
    #print("left {:.15f} {:.15f} {:.2f} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))
    #loc_right = Localization(grid=gri_right, dtm=None)
    #lonlatalt = loc_right.forward(row, col, h)
    #print("right {:.15f} {:.15f} {:.2f} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))



    #grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_left))
    #grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid.tif")



    matches_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches-crop.npy")
    matches = np.load(matches_filename)
    #print(matches.shape)
    #fname = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches-crop.tif")
    #with rio.open(fname, 'w', height=matches.shape[0], width=matches.shape[1], count=1, driver='GTiff', dtype=matches.dtype) as dst:
    #    dst.write_band(1, matches)
    #print("matches {}".format(matches[1:10,0:2]))

    point_ecef, point_wgs84, __ = epipolar_triangulation(matches, None,'sift',gri_left,gri_right,grid_left_filename,grid_right_filename)

    valid = [4584341.53073309455066919326782, 572313.952957414789125323295593, 4382784.34314652625471353530884]
    print(valid - point_ecef[0,:])
    assert(valid == pytest.approx(point_ecef[0,:],abs=1.0))




@pytest.mark.unit_tests
def test_epi_triangulation_sift_rpc():
    """
    Test epipolar triangulation
    """

    data_folder = test_path()
    id_scene = 'PHR1B_P_201709281038045_SEN_PRG_FC_178608-001'
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    geom_model_left = RPC.from_any(file_geom, topleftconvention=True)
    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    geom_model_right = RPC.from_any(file_geom, topleftconvention=True)

    #row = 10.0
    #col = 1000.0
    #h = 10.0
    #loc_left = Localization(grid=gri_left, dtm=None)
    #lonlatalt = loc_left.forward(row, col, h)
    #print("left {:.15f} {:.15f} {:.2f} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))
    #loc_right = Localization(grid=gri_right, dtm=None)
    #lonlatalt = loc_right.forward(row, col, h)
    #print("right {:.15f} {:.15f} {:.2f} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))



    #grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_left))
    #grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid.tif")



    matches_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches-crop.npy")
    matches = np.load(matches_filename)
    #print(matches.shape)
    #fname = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches-crop.tif")
    #with rio.open(fname, 'w', height=matches.shape[0], width=matches.shape[1], count=1, driver='GTiff', dtype=matches.dtype) as dst:
    #    dst.write_band(1, matches)

    point_ecef, point_wgs84, __ = epipolar_triangulation(matches, None,'sift',geom_model_left,geom_model_right,grid_left_filename,grid_right_filename)
    valid = [4584341.53073309455066919326782, 572313.952957414789125323295593, 4382784.34314652625471353530884]
    #print(valid - point_ecef[0,:])
    assert(valid == pytest.approx(point_ecef[0,:],abs=1e-3))



@pytest.mark.unit_tests
def test_epi_triangulation_sift_distance():
    """
    Test epipolar triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___,gri_right = prepare_loc('ellipsoide',id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc('ellipsoide',id_scene_left)

    #row = 10.0
    #col = 1000.0
    #h = 10.0
    #loc_left = Localization(grid=gri_left, dtm=None)
    #lonlatalt = loc_left.forward(row, col, h)
    #print("left {:.15f} {:.15f} {:.2f} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))
    #loc_right = Localization(grid=gri_right, dtm=None)
    #lonlatalt = loc_right.forward(row, col, h)
    #print("right {:.15f} {:.15f} {:.2f} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))



    #grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_left))
    #grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid.tif")



    matches_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches-crop.npy")
    matches = np.load(matches_filename)
    print(matches.shape)
    fname = os.path.join(os.environ["TESTPATH"], 'triangulation', "matches-crop.tif")
    #with rio.open(fname, 'w', height=matches.shape[0], width=matches.shape[1], count=1, driver='GTiff', dtype=matches.dtype) as dst:
    #    dst.write_band(1, matches)
    #print("matches {}".format(matches[1:10,0:2]))

    point_ecef, point_wgs84, residuals = epipolar_triangulation(matches, None,'sift',gri_left,gri_right,grid_left_filename,grid_right_filename, residues = True)

def stats_diff(cloud,array_epi):
    coords_x = cloud.x.values
    coords_y = cloud.y.values
    coords_z = cloud.z.values
    diff_x = abs(coords_x - array_epi[:,:,0])
    diff_y = abs(coords_y - array_epi[:, :, 1])
    diff_z = abs(coords_z - array_epi[:, :, 2])

    stats = np.zeros([3,3])
    mean_x = np.mean(diff_x)
    min_x = np.min(diff_x)
    max_x = np.max(diff_x)
    stats[0,:] = [mean_x,min_x,max_x]

    mean_y = np.mean(diff_y)
    min_y = np.min(diff_y)
    max_y = np.max(diff_y)
    stats[1,:] = [mean_y,min_y,max_y]
    mean_z = np.mean(diff_z)
    min_z= np.min(diff_z)
    max_z = np.max(diff_z)
    stats[2,:] = [mean_z,min_z,max_z]
    return stats



def create_dataset(disp,point_wgs84,point_ecef,residuals):
    array_shape = disp.disp.values.shape
    array_epi_wgs84 = point_wgs84.reshape((array_shape[0], array_shape[1], 3))
    array_epi_ecef = point_ecef.reshape((array_shape[0], array_shape[1], 3))
    array_residuals = residuals.reshape((array_shape[0], array_shape[1]))
    pc_dataset = xr.Dataset({'pc_wgs84_x': (['row', 'col'], array_epi_wgs84[:, :, 0]),
                             'pc_wgs84_y': (['row', 'col'], array_epi_wgs84[:, :, 1]),
                             'pc_wgs84_z': (['row', 'col'], array_epi_wgs84[:, :, 2]),
                             'pc_ecef_x': (['row', 'col'], array_epi_ecef[:, :, 0]),
                             'pc_ecef_y': (['row', 'col'], array_epi_ecef[:, :, 1]),
                             'pc_ecef_z': (['row', 'col'], array_epi_ecef[:, :, 2]),
                             'residues': (['row', 'col'], array_residuals)},
                            coords={'row': disp.coords['row'], 'col': disp.coords['col'] })
    return pc_dataset

@pytest.mark.unit_tests
def test_epi_triangulation_disp_rpc():
    """
     Test epipolar triangulation
    """
    data_folder = test_path()
    id_scene = 'PHR1B_P_201709281038045_SEN_PRG_FC_178608-001'
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    geom_model_left = RPC.from_any(file_geom, topleftconvention=True)
    id_scene = 'PHR1B_P_201709281038393_SEN_PRG_FC_178609-001'
    file_geom = os.path.join(data_folder, 'rpc/{}.geom'.format(id_scene))
    geom_model_right = RPC.from_any(file_geom, topleftconvention=True)



    #grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
    #                                  "grid_{}.tif".format(id_scene_left))
    #grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
    #                                   "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid.tif")

    disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)

    start = time.time()
    point_ecef, point_wgs84, residuals = epipolar_triangulation(disp, None, 'disp', geom_model_left, geom_model_right, grid_left_filename,
                                                   grid_right_filename, residues = True)
    end = time.time()
    print("elapsed time {}".format(end - start))
    pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)
    disp = xr.merge((disp, pc_dataset))
    out_disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "out_disparity_triangulation_rpc.nc")
    disp.to_netcdf(out_disp_filename)

    #open cloud
    cloud_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "cloud_ECEF.nc")
    cloud = xr.load_dataset(cloud_filename)
    array_shape = disp.disp.values.shape
    array_epi_ecef = point_ecef.reshape((array_shape[0], array_shape[1], 3))
    stats = stats_diff(cloud, array_epi_ecef)
    print(stats)
    index = 1492
    assert(point_ecef[index,0] == pytest.approx(cloud.x.values.flatten()[index],abs=1e-3))
    assert(point_ecef[index,1] == pytest.approx(cloud.y.values.flatten()[index],abs=1e-3))
    assert(point_ecef[index,2] == pytest.approx(cloud.z.values.flatten()[index],abs=1e-3))
    assert (stats[:,2] == pytest.approx([0,0,0], abs=6e-4))


@pytest.mark.unit_tests
def test_epi_triangulation_disp_rpc_roi():
    """
     Test epipolar triangulation
    """
    data_folder = test_path()
    file_geom = os.path.join(data_folder, 'rpc/phr_ventoux/left_image.geom')
    geom_model_left = RPC.from_any(file_geom, topleftconvention=True)
    file_geom = os.path.join(data_folder, 'rpc/phr_ventoux/right_image.geom')
    geom_model_right = RPC.from_any(file_geom, topleftconvention=True)


    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid_ventoux.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid_ventoux.tif")

    disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "disp1_ref.nc")
    disp = xr.load_dataset(disp_filename)

    start = time.time()
    point_ecef, point_wgs84, residuals = epipolar_triangulation(disp, None, 'disp', geom_model_left, geom_model_right, grid_left_filename,
                                                   grid_right_filename, residues = True)
    end = time.time()
    pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)

    #open cloud
    cloud_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "triangulation1_ref.nc")
    cloud = xr.load_dataset(cloud_filename)
    array_shape = disp.disp.values.shape
    array_epi_wgs84 = point_wgs84.reshape((array_shape[0], array_shape[1], 3))
    stats = stats_diff(cloud, array_epi_wgs84)
    #1492 first non masked index
    index = 100
    assert(point_wgs84[index,0] == pytest.approx(cloud.x.values.flatten()[index],abs=1e-8))
    assert(point_wgs84[index,1] == pytest.approx(cloud.y.values.flatten()[index],abs=1e-8))
    assert(point_wgs84[index,2] == pytest.approx(cloud.z.values.flatten()[index],abs=1e-3))
    assert (stats[:,2] == pytest.approx([0,0,0], abs=6e-4))




@pytest.mark.unit_tests
def test_epi_triangulation_disp_grid():
    """
     Test epipolar triangulation
    """
    id_scene_left = "P1BP--2017092838284574CP"
    id_scene_right = "P1BP--2017092838319324CP"
    ___, gri_right = prepare_loc('ellipsoide', id_scene_right)

    ___, gri_left = prepare_loc('ellipsoide', id_scene_left)

    #grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
    #                                  "grid_{}.tif".format(id_scene_left))
    #grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
    #                                   "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid.tif")

    disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)

    start = time.time()
    point_ecef, point_wgs84, residuals = epipolar_triangulation(disp, None, 'disp', gri_left, gri_right, grid_left_filename,
                                                   grid_right_filename, residues = True)
    end = time.time()
    print("elapsed time {}".format(end - start))
    pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)
    disp = xr.merge((disp, pc_dataset))
    out_disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "out_disparity_triangulation.nc")
    disp.to_netcdf(out_disp_filename)

    #open cloud
    cloud_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "cloud_ECEF.nc")
    cloud = xr.load_dataset(cloud_filename)
    array_shape = disp.disp.values.shape
    array_epi_ecef = point_ecef.reshape((array_shape[0], array_shape[1], 3))

    stats = stats_diff(cloud, array_epi_ecef)
    #1492 first non masked index
    index = 1492
    assert(point_ecef[index,0] == pytest.approx(cloud.x.values.flatten()[index],abs=1e-3))
    assert(point_ecef[index,1] == pytest.approx(cloud.y.values.flatten()[index],abs=0.6))
    assert(point_ecef[index,2] == pytest.approx(cloud.z.values.flatten()[index],abs=0.6))
    assert (stats[:,2] == pytest.approx([0,0,0], abs=0.6))




@pytest.mark.unit_tests
def test_epi_triangulation_disp_grid_masked():
    """
     Test epipolar triangulation
    """
    id_scene_left = "P1BP--2017092838284574CP"
    id_scene_right = "P1BP--2017092838319324CP"
    ___, gri_right = prepare_loc('ellipsoide', id_scene_right)

    ___, gri_left = prepare_loc('ellipsoide', id_scene_left)

    #grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
    #                                  "grid_{}.tif".format(id_scene_left))
    #grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids',
    #                                   "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(os.environ["TESTPATH"], 'rectification_grids', "right_epipolar_grid.tif")

    disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)
    mask_array = disp.msk.values
    start = time.time()
    point_ecef, point_wgs84, residuals = epipolar_triangulation(disp, mask_array, 'disp', gri_left, gri_right, grid_left_filename,
                                                   grid_right_filename, residues = True)
    end = time.time()
    print("elapsed time {}".format(end - start))
    #pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)
    #disp = xr.merge((disp, pc_dataset))
    #out_disp_filename = os.path.join(os.environ["TESTPATH"], 'triangulation', "out_disparity_triangulation_masked.nc")
    #disp.to_netcdf(out_disp_filename)
    assert (np.array_equal(point_ecef[0,:],[0,0,0]))
