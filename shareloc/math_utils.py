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

"""
This module contains the mathematical functions for shareloc
"""

#------------------------------------------------------------------------------
def interpol_bilin(mats,nl,nc,dl,dc):
        """
        bilinear interpolation on multi layer  matrix
        :param mats: multi layer grid (: , nl,nc)
        :param nl: line number of mats
        :param nc: column number of mats
        :param dl: position (line)
        :param dc: position (column)
        :return interpolated value on each layer
        :rtype list
        """
        if (dl < 0):
            i1 = 0
        elif (dl >= nl-1):
            i1 = nl - 2
        else:
            i1 = int(np.floor(dl))
        i2 = i1+1

        if (dc < 0):
            j1 = 0
        elif (dc >= nc-1):
            j1 = nc - 2
        else:
            j1 = int(np.floor(dc))
        j2 = j1+1
        #(u,v) are subpixel distance to interpolate along each axis
        u = dc - j1
        v = dl - i1
        #Altitude
        matis=[]
        for mat in mats:
            mati = (1-u)*(1-v)*mat[:,i1,j1] + u*(1-v)*mat[:,i1,j2] +\
                   (1-u)*v*mat[:,i2,j1]     + u*v*mat[:,i2,j2]
            matis.append(mati)
        return matis
