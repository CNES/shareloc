#!/usr/bin/env python
# coding: utf8
#
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
#
"""
Test module for triangulation class shareloc/geofunctions/triangulation.py
"""

# Standard imports
import os

# Third party imports
import numpy as np
import pytest
import xarray as xr

# Shareloc imports
from shareloc.geofunctions.triangulation import distance_point_los, epipolar_triangulation, sensor_triangulation
from shareloc.geomodels.grid import Grid
from shareloc.geomodels.rpc import RPC

# Shareloc test imports
from ..helpers import data_path


def prepare_loc(alti="geoide", id_scene="P1BP--2017030824934340CP"):
    """
    Read multiH grid

    :param alti: alti validation dir
    :param id_scene: scene ID
    :return: multi H grid
    :rtype: str
    """
    data_folder = data_path(alti, id_scene)
    # load grid
    gld = os.path.join(data_folder, f"GRID_{id_scene}.tif")
    gri = Grid(gld)

    return gri


@pytest.mark.parametrize("col,row,h", [(1000.5, 1500.5, 10.0)])
@pytest.mark.unit_tests
def test_sensor_triangulation(row, col, h):
    """
    Test sensor triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    gri_right = prepare_loc("ellipsoide", id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    gri_left = prepare_loc("ellipsoide", id_scene_left)
    # init des predicteurs
    gri_right.estimate_inverse_loc_predictor()
    lonlatalt = gri_left.direct_loc_h(row, col, h)

    inv_row, inv_col, __ = gri_right.inverse_loc(lonlatalt[0], lonlatalt[1], lonlatalt[2])

    matches = np.zeros([1, 4])
    matches[0, :] = [col, row, inv_col, inv_row]
    # matches[1,:] = [lig + 10, col + 5, inv_lig + 12, inv_col + 7]

    __, point_wgs84, distance = sensor_triangulation(matches, gri_left, gri_right, residues=True)

    assert lonlatalt[0] == pytest.approx(point_wgs84[0, 0], abs=1e-8)
    assert lonlatalt[1] == pytest.approx(point_wgs84[0, 1], abs=1e-8)
    assert lonlatalt[2] == pytest.approx(point_wgs84[0, 2], abs=8e-3)
    assert distance == pytest.approx(0.0, abs=1e-3)


@pytest.mark.unit_tests
def test_triangulation_residues():
    """
    Test triangulation residues on simulated LOS
    """

    class SimulatedLOS:
        """line of sight class"""

        def __init__(self):
            self.sis = np.array([[100.0, 10.0, 200.0], [100.0, 10.0, 200.0]])
            self.vis = np.array([[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]])

        def print_sis(self):
            """
            print los hat
            """
            print(self.sis)

        def print_vis(self):
            """
            print los viewing vector
            """
            print(self.vis)

    los = SimulatedLOS()

    distance = 10.0
    point = los.sis + 100.0 * los.vis + distance * np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
    residue = distance_point_los(los, point)
    assert distance == pytest.approx(residue, abs=1e-9)


@pytest.mark.unit_tests
def test_epi_triangulation_sift():
    """
    Test epipolar triangulation
    """
    id_scene_right = "P1BP--2017092838319324CP"
    gri_right = prepare_loc("ellipsoide", id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    gri_left = prepare_loc("ellipsoide", id_scene_left)

    grid_left_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(data_path(), "rectification_grids", "right_epipolar_grid.tif")

    matches_filename = os.path.join(data_path(), "triangulation", "matches-crop.npy")
    matches = np.load(matches_filename)

    point_ecef, __, __ = epipolar_triangulation(
        matches, None, "sift", gri_left, gri_right, grid_left_filename, grid_right_filename
    )
    valid = [4584341.37359843123704195022583, 572313.675204274943098425865173, 4382784.51356450468301773071289]
    assert valid == pytest.approx(point_ecef[0, :], abs=0.5)


@pytest.mark.unit_tests
def test_epi_triangulation_sift_rpc():
    """
    Test epipolar triangulation
    """

    data_folder = data_path()
    id_scene = "PHR1B_P_201709281038045_SEN_PRG_FC_178608-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    geom_model_left = RPC.from_any(file_geom, topleftconvention=True)
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    print(file_geom)
    geom_model_right = RPC.from_any(file_geom, topleftconvention=True)

    grid_left_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(data_path(), "rectification_grids", "right_epipolar_grid.tif")

    matches_filename = os.path.join(data_path(), "triangulation", "matches-crop.npy")
    matches = np.load(matches_filename)

    point_ecef, __, __ = epipolar_triangulation(
        matches, None, "sift", geom_model_left, geom_model_right, grid_left_filename, grid_right_filename
    )
    valid = [4584341.37359843123704195022583, 572313.675204274943098425865173, 4382784.51356450468301773071289]
    # print(valid - point_ecef[0,:])
    assert valid == pytest.approx(point_ecef[0, :], abs=1e-3)


def stats_diff(cloud, array_epi):
    """
    compute difference statistics between dataset and shareloc results
    :param cloud :  CARS dataset
    :type cloud : xarray dataset
    :param array_epi :  shareloc array
    :type array_epi : numpy.array
    :return stats [mean,min, max]
    :rtype numpy.array
    """
    coords_x = cloud.x.values
    coords_y = cloud.y.values
    coords_z = cloud.z.values
    diff_x = abs(coords_x - array_epi[:, :, 0])
    diff_y = abs(coords_y - array_epi[:, :, 1])
    diff_z = abs(coords_z - array_epi[:, :, 2])

    stats = np.zeros([3, 3])
    mean_x = np.mean(diff_x)
    min_x = np.min(diff_x)
    max_x = np.max(diff_x)
    stats[0, :] = [mean_x, min_x, max_x]

    mean_y = np.mean(diff_y)
    min_y = np.min(diff_y)
    max_y = np.max(diff_y)
    stats[1, :] = [mean_y, min_y, max_y]
    mean_z = np.mean(diff_z)
    min_z = np.min(diff_z)
    max_z = np.max(diff_z)
    stats[2, :] = [mean_z, min_z, max_z]
    return stats


def create_dataset(disp, point_wgs84, point_ecef, residuals):
    """
    create new dataset from existing one
    :param point_wgs84 :  points WGS84
    :type point_wgs84 : numpy.array
    :param point_ecef :  points ECEF
    :type point_ecef : numpy.array
    :param residuals :  traingulations residuals
    :type residuals : numpy.array
    :return dataset
    :rtype xarray.dataset
    """
    array_shape = disp.disp.values.shape
    array_epi_wgs84 = point_wgs84.reshape((array_shape[0], array_shape[1], 3))
    array_epi_ecef = point_ecef.reshape((array_shape[0], array_shape[1], 3))
    array_residuals = residuals.reshape((array_shape[0], array_shape[1]))
    pc_dataset = xr.Dataset(
        {
            "pc_wgs84_x": (["row", "col"], array_epi_wgs84[:, :, 0]),
            "pc_wgs84_y": (["row", "col"], array_epi_wgs84[:, :, 1]),
            "pc_wgs84_z": (["row", "col"], array_epi_wgs84[:, :, 2]),
            "pc_ecef_x": (["row", "col"], array_epi_ecef[:, :, 0]),
            "pc_ecef_y": (["row", "col"], array_epi_ecef[:, :, 1]),
            "pc_ecef_z": (["row", "col"], array_epi_ecef[:, :, 2]),
            "residues": (["row", "col"], array_residuals),
        },
        coords={"row": disp.coords["row"], "col": disp.coords["col"]},
    )
    return pc_dataset


@pytest.mark.unit_tests
def test_epi_triangulation_disp_rpc():
    """
    Test epipolar triangulation
    """
    data_folder = data_path()
    id_scene = "PHR1B_P_201709281038045_SEN_PRG_FC_178608-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    geom_model_left = RPC.from_any(file_geom, topleftconvention=True)
    id_scene = "PHR1B_P_201709281038393_SEN_PRG_FC_178609-001"
    file_geom = os.path.join(data_folder, f"rpc/{id_scene}.geom")
    geom_model_right = RPC.from_any(file_geom, topleftconvention=True)

    # grid_left_filename = os.path.join(data_path(), "rectification_grids",
    #                                  "grid_{}.tif".format(id_scene_left))
    # grid_right_filename = os.path.join(data_path(), "rectification_grids",
    #                                   "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(data_path(), "rectification_grids", "right_epipolar_grid.tif")

    disp_filename = os.path.join(data_path(), "triangulation", "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)

    point_ecef, point_wgs84, residuals = epipolar_triangulation(
        disp, None, "disp", geom_model_left, geom_model_right, grid_left_filename, grid_right_filename, residues=True
    )
    pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)
    disp = xr.merge((disp, pc_dataset))

    # open cloud
    cloud_filename = os.path.join(data_path(), "triangulation", "cloud_ECEF.nc")
    cloud = xr.load_dataset(cloud_filename)
    array_shape = disp.disp.values.shape
    array_epi_ecef = point_ecef.reshape((array_shape[0], array_shape[1], 3))
    stats = stats_diff(cloud, array_epi_ecef)
    index = 1492
    assert point_ecef[index, 0] == pytest.approx(cloud.x.values.flatten()[index], abs=1e-3)
    assert point_ecef[index, 1] == pytest.approx(cloud.y.values.flatten()[index], abs=1e-3)
    assert point_ecef[index, 2] == pytest.approx(cloud.z.values.flatten()[index], abs=1e-3)
    assert stats[:, 2] == pytest.approx([0, 0, 0], abs=6e-4)


@pytest.mark.unit_tests
def test_epi_triangulation_disp_rpc_roi():
    """
    Test epipolar triangulation
    """
    data_folder = data_path()
    file_geom = os.path.join(data_folder, "rpc/phr_ventoux/left_image.geom")
    geom_model_left = RPC.from_any(file_geom, topleftconvention=True)
    file_geom = os.path.join(data_folder, "rpc/phr_ventoux/right_image.geom")
    geom_model_right = RPC.from_any(file_geom, topleftconvention=True)

    grid_left_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid_ventoux.tif")
    grid_right_filename = os.path.join(data_path(), "rectification_grids", "right_epipolar_grid_ventoux.tif")

    disp_filename = os.path.join(data_path(), "triangulation", "disp1_ref.nc")
    disp = xr.load_dataset(disp_filename)

    __, point_wgs84, __ = epipolar_triangulation(
        disp,
        None,
        "disp",
        geom_model_left,
        geom_model_right,
        grid_left_filename,
        grid_right_filename,
        residues=True,
        fill_nan=True,
    )
    # pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)

    # open cloud
    cloud_filename = os.path.join(data_path(), "triangulation", "triangulation1_ref.nc")
    cloud = xr.load_dataset(cloud_filename)
    array_shape = disp.disp.values.shape
    array_epi_wgs84 = point_wgs84.reshape((array_shape[0], array_shape[1], 3))
    stats = stats_diff(cloud, array_epi_wgs84)
    # 1492 first non masked index
    index = 100
    assert point_wgs84[index, 0] == pytest.approx(cloud.x.values.flatten()[index], abs=1e-8)
    assert point_wgs84[index, 1] == pytest.approx(cloud.y.values.flatten()[index], abs=1e-8)
    assert point_wgs84[index, 2] == pytest.approx(cloud.z.values.flatten()[index], abs=1e-3)
    assert stats[:, 2] == pytest.approx([0, 0, 0], abs=6e-4)


@pytest.mark.unit_tests
def test_epi_triangulation_disp_grid():
    """
    Test epipolar triangulation
    """
    id_scene_left = "P1BP--2017092838284574CP"
    id_scene_right = "P1BP--2017092838319324CP"
    gri_right = prepare_loc("ellipsoide", id_scene_right)

    gri_left = prepare_loc("ellipsoide", id_scene_left)

    # grid_left_filename = os.path.join(data_path(), "rectification_grids",
    #                                  "grid_{}.tif".format(id_scene_left))
    # grid_right_filename = os.path.join(data_path(), "rectification_grids",
    #                                   "grid_{}.tif".format(id_scene_right))

    grid_left_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(data_path(), "rectification_grids", "right_epipolar_grid.tif")

    disp_filename = os.path.join(data_path(), "triangulation", "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)

    point_ecef, point_wgs84, residuals = epipolar_triangulation(
        disp, None, "disp", gri_left, gri_right, grid_left_filename, grid_right_filename, residues=True
    )
    pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)
    disp = xr.merge((disp, pc_dataset))

    # open cloud
    cloud_filename = os.path.join(data_path(), "triangulation", "cloud_ECEF.nc")
    cloud = xr.load_dataset(cloud_filename)
    array_shape = disp.disp.values.shape
    array_epi_ecef = point_ecef.reshape((array_shape[0], array_shape[1], 3))

    stats = stats_diff(cloud, array_epi_ecef)
    # 1492 first non masked index
    index = 1492
    assert point_ecef[index, 0] == pytest.approx(cloud.x.values.flatten()[index], abs=0.3)
    assert point_ecef[index, 1] == pytest.approx(cloud.y.values.flatten()[index], abs=0.3)
    assert point_ecef[index, 2] == pytest.approx(cloud.z.values.flatten()[index], abs=0.3)
    assert stats[:, 2] == pytest.approx([0, 0, 0], abs=0.5)


@pytest.mark.unit_tests
def test_epi_triangulation_disp_grid_masked():
    """
    Test epipolar triangulation
    """
    id_scene_left = "P1BP--2017092838284574CP"
    id_scene_right = "P1BP--2017092838319324CP"
    gri_right = prepare_loc("ellipsoide", id_scene_right)

    gri_left = prepare_loc("ellipsoide", id_scene_left)

    grid_left_filename = os.path.join(data_path(), "rectification_grids", "left_epipolar_grid.tif")
    grid_right_filename = os.path.join(data_path(), "rectification_grids", "right_epipolar_grid.tif")

    disp_filename = os.path.join(data_path(), "triangulation", "disparity-crop.nc")
    disp = xr.load_dataset(disp_filename)
    mask_array = disp.msk.values
    point_ecef, __, __ = epipolar_triangulation(
        disp, mask_array, "disp", gri_left, gri_right, grid_left_filename, grid_right_filename, residues=True
    )
    # pc_dataset = create_dataset(disp, point_wgs84, point_ecef, residuals)
    # disp = xr.merge((disp, pc_dataset))
    assert np.array_equal(point_ecef[0, :], [0, 0, 0])
