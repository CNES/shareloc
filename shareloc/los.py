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
from shareloc.proj_utils import coordinates_conversion

class LOS:
    """ line of sight class
    """
    def __init__(self, sensor_positions, geometrical_model, alt_min_max  = None):
        """
        Constructor
        :param sensor_positions: sensor_positions
        :type sensor_positions: numpy array (Nx2)
        :param geometrical_model: geometrical model
        :type geometrical_model: shareloc.grid or shareloc..rpc.rpc
        :param alt_min_max: min/max  altitude to compute los, if None model min/max will be used
        :type alt_min_max : list

        """

        self.geometrical_model = geometrical_model
        self.sensors_positions = sensor_positions
        self.los_creation(alt_min_max)

    def los_creation(self,alt_min_max):
        """
        create los from extrema : los starting point, and normalized viewing vector
        :param alt_min_max: min/max  altitude to compute los, if None model min/max will be used
        :type alt_min_max : list
        """

        self.los_nb = self.sensors_positions.shape[0]
        if(alt_min_max is None):
            (alt_min,alt_max) = self.geometrical_model.get_alt_min_max()

        # los construction right
        los_extrema = np.zeros([2 * self.los_nb, 3])
        list_col, list_row = (self.sensors_positions[:, 0], self.sensors_positions[:, 1])
        los_extrema[np.arange(0, 2 *self.los_nb, 2), :] = self.geometrical_model.direct_loc_h(list_row, list_col, alt_max)
        los_extrema[np.arange(1, 2 * self.los_nb, 2), :] = self.geometrical_model.direct_loc_h(list_row, list_col, alt_min)

        in_crs = 4326
        out_crs = 4978
        ecef_coord = coordinates_conversion(los_extrema, in_crs, out_crs)
        self.sis = ecef_coord[0::2, :]
        vis = self.sis - ecef_coord[1::2, :]
        # /!\ normalized
        #
        #  direction vector creation
        vis_norm = np.linalg.norm(vis, axis=1)
        rep_vis_norm = np.tile(vis_norm,(3,1)).transpose()
        self.vis = vis / rep_vis_norm



