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

import numpy as np
from shareloc.rectification.rectification_grid import rectification_grid
from shareloc.proj_utils import coordinates_conversion
from shareloc.los import LOS

def sensor_triangulation(matches, geometrical_model_left,geometrical_model_right,left_min_max = None,
                         right_min_max = None, residues = False, fill_nan = False):
    """
    triangulation in sensor geometry

    according to the formula:
    .. math::
        x =
        \\left(\\sum_i I-\\hat v_i \\hat v_i^\\top\\right)^{-1} \\left(\\sum_i (I-\\hat v_i \\hat v_i^\\top) s_i\\right)
    Delvit J.M. et al. "The geometric supersite of Salon de Provence", ISPRS Congress Paris, 2006.


    :param matches :  matches in sensor coordinates Nx[row (left), col (left), row (right), col (right)]
    :type matches : np.array
    :param geometrical_model_left : left image geometrical model
    :type geometrical_model_left : shareloc.grid or shareloc.rpc
    :param geometrical_model_right : right image geometrical model
    :type geometrical_model_right : shareloc.grid or shareloc.rpc
    :param left_min_max : left min/max for los creation, if None model min/max will be used
    :type left_min_max : list
    :param right_min_max : right min/max for los creation, if None model min/max will be used
    :type right_min_max : list
    :param residues : calculates residues (distance in meters between los and 3D points)
    :type residues : boolean
    :param fill_nan : fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
        otherwise
    :type fill_nan : boolean
    :return intersections in cartesian crs, intersections in wgs84 crs and optionnaly residues
    :rtype (numpy.array,numpy,array,numpy.array)
    """
    #los construction
    matches_left = matches[:,0:2]
    left_los = LOS(matches_left, geometrical_model_left, left_min_max, fill_nan)
    matches_right = matches[:, 2:4]
    right_los = LOS(matches_right, geometrical_model_right, right_min_max, fill_nan)

    #los conversion
    #los intersection
    intersections_ecef = los_triangulation(left_los,right_los)
    in_crs = 4978
    out_crs = 4326
    intersections_wgs84 = coordinates_conversion(intersections_ecef, in_crs, out_crs)
    #refine matches
    intersections_residues = None
    if residues is True:
        intersections_residues = distance_point_los(left_los,intersections_ecef)
    return intersections_ecef,intersections_wgs84, intersections_residues


def distance_point_los(los, points):
    """
    distance between points and los
    norm of cross prodcut between vector defined by los sis and point and los vis

    :param los :  line of sight
    :type los : shareloc.los
    :param points :  3D points
    :type points : numpy.array
    :return distance
    :rtype numpy.array
    """
    # norm(BA vect u)/norm(u)

    vis = los.vis
    vect_sis_p = points - los.sis
    dist = np.linalg.norm(np.cross(vect_sis_p, vis), axis = 1)
    return dist


def los_triangulation(left_los,right_los):
    """
    los triangulation

    :param left_los :  left los
    :type left_los : shareloc.los
    :param right_los :  right los
    :type right_los : shareloc.los
    :return intersections in cartesian crs
    :rtype numpy.array
    """
    vis = np.dstack((left_los.vis, right_los.vis))
    vis = np.swapaxes(vis,1,2)


    sis = np.dstack((left_los.sis, right_los.sis))
    sis = np.swapaxes(sis, 1, 2)

    vivi = vis[..., :, np.newaxis] * vis[..., np.newaxis, :]
    id_vivi = np.eye(3) - vivi
    sum_id_vivi = np.nansum(id_vivi, axis=-3)

    id_vivi_si = np.nansum(id_vivi * sis[..., np.newaxis, :], axis=-1)
    sum_id_vivi_si = np.nansum(id_vivi_si, axis=-2)

    inv_sum_id_vivi = np.linalg.inv(sum_id_vivi)
    intersection = np.nansum(inv_sum_id_vivi * sum_id_vivi_si[..., np.newaxis, :], axis=-1)
    return intersection


def transform_disp_to_matches(disp, mask = None):
    """
    transform disparity map to matches

    :param disp :  dispartiy xarray
    :type disp : xarray
    :param mask :  mask
    :type mask : numpy.array
    :return epipolar matches and logical non masked values
    :rtype list of numpy.arrray
    """


    col,row  = np.meshgrid(disp.col,disp.row)
    disp_array = disp.disp.values
    if mask is not None:
        values_ok = mask > 0
        col = col[values_ok]
        row = row[values_ok]
        disp_array = disp_array[values_ok]
    else:
        values_ok = np.full_like(disp_array, True, dtype=bool)
    epi_left_pos  = np.vstack((col.flatten(),row.flatten()))

    epi_right_pos = np.vstack((col.flatten() + disp_array.flatten(),row.flatten()))
    return (epi_left_pos.transpose(),epi_right_pos.transpose(),values_ok.flatten())




def epipolar_triangulation(matches,mask, matches_type, geometrical_model_left,geometrical_model_right,
                           grid_left,grid_right,left_min_max = None,right_min_max = None,residues = False,
                           fill_nan = False):
    """
    epipolar triangulation

    :param matches :  matches
    :type matches :
    :param mask :  mask
    :type mask :
    :param matches_type :  'disp' or 'sift'
    :type matches_type : str
    :param geometrical_model_left : left image geometrical model
    :type geometrical_model_left : shareloc.grid or shareloc.rpc
    :param geometrical_model_right : right image geometrical model
    :type geometrical_model_right : shareloc.grid or shareloc.rpc
    :param grid_left : left rectification grid filename
    :type grid_left : str
    :param grid_right : right rectification grid filename
    :type grid_right : str
    :param left_min_max : left min/max for los creation, if None model min/max will be used
    :type left_min_max : list
    :param right_min_max : right min/max for los creation, if None model min/max will be used
    :type right_min_max : list
    :param residues : calculates residues (distance in meters)
    :type residues : boolean
    :param fill_nan : fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
        otherwise
    :type fill_nan : boolean
    :return intersections in cartesian crs
    :rtype numpy.array
    """

    # retrieve point matches in sensor geometry
    if matches_type is 'sift':
        epi_pos_left = matches[:,0:2]
        epi_pos_right= matches[:,2:4]
        values_ok = np.full((matches.shape[0]), True, dtype=bool)
    elif matches_type is 'disp':
        [epi_pos_left, epi_pos_right, values_ok] = transform_disp_to_matches(matches,mask)
        #epi_pos_left = matches[:,0:2]
        #epi_disp_right= matches[:,2:4]

    else:
            raise Exception(
                'matches type should be sift or disp')

    #interpolate left
    rectif_grid_left = rectification_grid(grid_left)
    rectif_grid_right = rectification_grid(grid_right)
    #interpolate_right
    #print("epi  left 0 {}".format(epi_pos_left[0,:]))
    #print("epi right 0 {}".format(epi_pos_right[0,:]))
    matches_sensor_left = rectif_grid_left.interpolate(epi_pos_left)
    matches_sensor_right = rectif_grid_right\
        .interpolate(epi_pos_right)
    #print("epi  left 0 {} sensor 0 {}".format(epi_pos_left[0,:], matches_sensor_left[0,:]))
    #print("epi right 0 {} sensor 0 {}".format(epi_pos_right[0,:], matches_sensor_right[0,:]))
    matches_sensor = np.concatenate((matches_sensor_left,matches_sensor_right),axis = 1)
    #print("matches shape {}".format(matches_sensor.shape))
    #print("matches  {}".format(matches_sensor[0,:]))
    # triangulate positions in sensor geometry
    [intersections_ecef, intersections_wgs84, intersections_residues] = sensor_triangulation(matches_sensor,
                                        geometrical_model_left, geometrical_model_right,
                                        left_min_max, right_min_max, residues, fill_nan)

    tab_size = values_ok.shape[0]
    intersections_ecef_masked = np.zeros((tab_size,3))
    intersections_wgs84_masked = np.zeros((tab_size, 3))
    intersections_residues_masked = np.zeros((tab_size, 1))
    intersections_ecef_masked[values_ok,:] = intersections_ecef
    intersections_wgs84_masked[values_ok, :] = intersections_wgs84
    intersections_residues_masked[values_ok,0] = intersections_residues
    return intersections_ecef_masked, intersections_wgs84_masked, intersections_residues_masked
