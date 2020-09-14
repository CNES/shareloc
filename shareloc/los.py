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


class LOS:
    """ line of sight class
    """
    def __init__(self, sensor_positions, geometrical_model, alt_min = 0.0 , alt_max = 4000.0):
        """
        Constructor
        :param sensor_positions: sensor_positions
        :type sensor_positions: numpy array (Nx2)
        :param geometrical_model: geometrical model
        :type geometrical_model: shareloc.grid or shareloc..rpc.rpc
        :param alt_min: min altitude to compute los
        :type alt_min: float

        """
        self.alt_min = alt_min
        self.alt_max = alt_max
        self.geometrical_model = geometrical_model
        self.sensors_positions = sensor_positions
        self.sis, self.vis = geometrical_model.los_creation(sensor_positions, alt_min, alt_max)

    def intersection(sis, vis, get_index_valid=False):
        """
        Find intersection
        Args:
            sis: lines of sights starting points
            vis: lines of sights direction vectors /!\ normalized

            sis[m, i, [x, y, z]] and vis[m, i, [x, y, z]]
            with m: the mth match
                 i: the ith corresponding image/model
                 [x, y, z]: coordinates of the starting point (ecef) for sis
                            coordinates if the normalized direction vector
            get_index_valid: index of valid matches (number of lines > 1)
        Returns:
            intersection: position of the closest points (numpy array)
            intersection[m, [x, y, z]]
            with m: the mth match
                 [x, y, z]: coordinates of the intersection
            index_valid if get_index_valid is True

        according to the formula:
        x= \left(\sum_i I-\hat v_i \hat v_i^\top\right)^{-1} \left(\sum_i (I-\hat v_i \hat v_i^\top) s_i\right)
        Delvit J.M. et al. "The geometric supersite of Salon de Provence", ISPRS Congress Paris, 2006.

        """
        # remove unique lines
        nb_los = np.ravel(np.nansum(np.isfinite(sis[..., 0]), axis=-1))
        index_valid = np.where(nb_los > 1)
        sis_no_nan = sis[index_valid]
        vis_no_nan = vis[index_valid]

        # compute intersection
        vivi = vis_no_nan[..., :, np.newaxis] * vis_no_nan[..., np.newaxis, :]
        id_vivi = np.eye(3) - vivi
        sum_id_vivi = np.nansum(id_vivi, axis=-3)

        id_vivi_si = np.nansum(id_vivi * sis_no_nan[..., np.newaxis, :], axis=-1)
        sum_id_vivi_si = np.nansum(id_vivi_si, axis=-2)

        inv_sum_id_vivi = np.linalg.inv(sum_id_vivi)
        intersection = np.nansum(inv_sum_id_vivi * sum_id_vivi_si[..., np.newaxis, :], axis=-1)

        if get_index_valid is True:
            return intersection, index_valid
        else:
            return intersection

