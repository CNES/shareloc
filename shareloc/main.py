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

import os.path as op
from shareloc.gld_xh import fct_coloc, gld_xH
from shareloc.mnt import mnt
import os


def test_path():
    """
    return the data folder
    :return: data path.
    :rtype: str
    """    
    data_folder = op.join(os.environ["TESTPATH"])
    return data_folder




if __name__ == '__main__':
    """
    calcul_gld3d -all -path_visee ./grilles_gld_xH -nom_visee P1BP--2017030824934340CP_H1 -type_visee LocalisateurGrille_Directe \
    -mnt ./MNT_extrait/mnt_extrait -repter GRS80:G-D/:H-M -path_grille . -nom_grille test_intersect_euclide -convention BABEL \
    -format BSQ -nbcol 200 -nblig 200 -pascol 50 -paslig 60 -j0 0 -i0 0 -col0 100.5 -lig0 100.5 -matrice 1 0 0 1
    """

    import time
    data_folder = test_path()
    
    #chargement du mnt
    fic = op.join(data_folder,'MNT_extrait/mnt_extrait.c1')
    mntbsq = mnt(fic)
    
    #chargement des grilles
    gld = op.join(data_folder,'grilles_gld_xH/P1BP--2017030824934340CP_H1.hd')
    gri = gld_xH(gld)
    
    #init des predicteurs
    gri.init_pred_loc_inv()
    #
    start_time = time.time()
    gricol = fct_coloc(gri, gri, mntbsq,0.5, 0.5, 10, 100, 20, 20)
    gri_gld = gri.fct_gld_mnt(100.5, 100.5, 60, 50, 200, 200, mntbsq)
    interval = time.time() - start_time
    print('Total time in seconds 200x200 coloc: {:.2f}s'.format(interval))
    
    
    #validation lecture de grille

