# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 19:07:16 2020

@author: gresloud
"""

import os
import pytest
from utils import test_path

from shareloc.rpc.rpc_phr_v2 import FonctRatD


def test_rpc_eucl():
    data_folder = test_path()

    fichier_direct_euclide = os.path.join(data_folder,'rpc/rpc_dir.euc')
    fichier_inverse_euclide =  os.path.join(data_folder,'rpc/rpc_inv.euc')

    fctrat = FonctRatD(fichier_direct_euclide, fichier_inverse_euclide)

    col = 600
    lig = 200
    alt = 125

    (lon,lat) = fctrat.evalue_loc_d(col,lig,alt)
    print (lon,lat)

    (col_ar,lig_ar) = fctrat.evalue_loc_i(lon,lat,alt)
    print (col_ar,lig_ar)
    assert(col_ar == pytest.approx(col,4))
    assert(lig_ar == pytest.approx(lig,4))

def test_rpc_phrdimap():
    data_folder = test_path()

    fichier_dimap = os.path.join(data_folder,'rpc/PHRDIMAP.XML')

    fctrat = FonctRatD(fichier_dimap)

    col = 600
    lig = 200
    alt = 125

    (lon,lat) = fctrat.evalue_loc_d(col,lig,alt)
    print (lon,lat)

    (col_ar,lig_ar) = fctrat.evalue_loc_i(lon,lat,alt)
    print (col_ar,lig_ar)
    assert(col_ar == pytest.approx(col,4))
    assert(lig_ar == pytest.approx(lig,4))
