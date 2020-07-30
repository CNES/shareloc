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

#-------------------------------------------------------------------------------




class sensor:
    """ sensor class
    """
    def __init__(self, mnt = None, grid = None, rpc = None ):
        self.mnt = mnt
        self.grid = grid
        self.rpc = rpc
        self.use_rpc = False
        if self.rpc is not None and self.grid is None:
            self.use_rpc = True


    def forward(self, lig, col, h):
        if self.use_rpc == True:
            print('use rpc')
            return self.rpc.evalue_loc_d(col,lig, h)
        else:
            return self.grid.fct_locdir_h(lig,col,h)

    def forward_dtm(self, lig, col):
        if self.use_rpc == True:
            print('forward dtm not yet impelemented for RPC model')
            return None
        else:
            if self.mnt is not None:
                return self.grid.fct_locdir_mnt(lig,col, self.mnt)
            else:
                print('forward dtm needs a mntl')
                return None

    def inverse(self,lat,lon,h):
        if self.use_rpc == True:
            print('use rpc')
            return self.rpc.evalue_loc_i(lon,lat,h)
        else:
            return self.grid.fct_locinv([lon,lat,h])



