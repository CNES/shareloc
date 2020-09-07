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
def lit_header_hdbabel(fic_mnt_bsq):
    path_header = op.dirname(fic_mnt_bsq)
    name_header = op.basename(fic_mnt_bsq).split('.')[0]+'.hd_babel'
    entete = {}
    hd_babel = {}
    with open(op.join(path_header,name_header),'r') as f:
        txt_header = f.readlines()
    for lig_header in txt_header:
        if lig_header.startswith('>>'):
            ligsplit = lig_header.split('\t')
            entete[ligsplit[1]] = ligsplit[2]

    hd_babel['x0'] = float(entete['LON_REF'])
    hd_babel['y0'] = float(entete['LAT_REF'])
    hd_babel['nc'] = int(entete['NB_LON'])
    hd_babel['nl'] = int(entete['NB_LAT'])
    hd_babel['px'] = float(entete['PAS_LON'])
    hd_babel['py'] = float(entete['PAS_LAT'])
    if int(entete['TYPE_CODE'])==2:
        hd_babel['codage'] = 'int16'
    return hd_babel

#------------------------------------------------------------------------------
def lit_grille_bsq(fic_bsq,nl,nc,codage):
    grille = np.fromfile(fic_bsq,dtype=codage).reshape((nl,nc))
    return grille
#------------------------------------------------------------------------------
def lit_hd_bsq(fic_hd, balises):
    """balises est un dico
    dico_out['nom'] = ('nom_bsq',type)
    """
    dico_out = {}
    with open(op.join(fic_hd),'r') as f:
        txt_header = f.readlines()
    entete = {}

    for i in range(0,txt_header.__len__(),2):
        entete[txt_header[i].strip()] = txt_header[i+1].strip()

    for var in balises:
        (nom, form) = balises[var]
        dico_out[var] = form(entete[nom])

    return dico_out
