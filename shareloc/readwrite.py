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
import os.path as op

"""
This module contains the input/output for shareloc
"""


#------------------------------------------------------------------------------
def read_hdbabel_header(fic_mnt_bsq):
    path_header = op.dirname(fic_mnt_bsq)
    name_header = op.basename(fic_mnt_bsq).split('.')[0]+'.hd_babel'
    header = {}
    hd_babel = {}
    with open(op.join(path_header,name_header),'r') as f:
        txt_header = f.readlines()
    for lig_header in txt_header:
        if lig_header.startswith('>>'):
            ligsplit = lig_header.split('\t')
            header[ligsplit[1]] = ligsplit[2]

    hd_babel['x0'] = float(header['LON_REF'])
    hd_babel['y0'] = float(header['LAT_REF'])
    hd_babel['nc'] = int(header['NB_LON'])
    hd_babel['nl'] = int(header['NB_LAT'])
    hd_babel['px'] = float(header['PAS_LON'])
    hd_babel['py'] = float(header['PAS_LAT'])
    if int(header['TYPE_CODE'])==2:
        hd_babel['data_type'] = 'int16'
    return hd_babel

#------------------------------------------------------------------------------
def read_bsq_grid(fic_bsq, nl, nc, data_type):
    grid = np.fromfile(fic_bsq,dtype=data_type).reshape((nl,nc))
    return grid
#------------------------------------------------------------------------------
def read_bsq_hd(fic_hd, tag):
    """balises est un dico
    dico_out['nom'] = ('nom_bsq',type)
    """
    dico_out = {}
    with open(op.join(fic_hd),'r') as f:
        txt_header = f.readlines()
    header = {}

    for i in range(0,txt_header.__len__(),2):
        header[txt_header[i].strip()] = txt_header[i+1].strip()

    for var in tag:
        (nom, form) = tag[var]
        dico_out[var] = form(header[nom])

    return dico_out
