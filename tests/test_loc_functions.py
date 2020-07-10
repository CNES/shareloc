#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright: (c) 2020 Centre National d'Etudes Spatiales

import os
import pytest
import shareloc.pyi3D_v3 as loc


def test_path():
    """
    return the data folder
    :return: data path.
    :rtype: str
    """    
    data_folder = os.path.join(os.environ["TESTPATH"])
    return data_folder

def prepare_loc():
    """
    Read multiH grid
    :return: multi H grid
    :rtype: str
    """   
    data_folder = test_path()
    
    #chargement du mnt
    fic = os.path.join(data_folder,'MNT_extrait/mnt_extrait.c1')
    mntbsq = loc.mnt(fic)
    
    #chargement des grilles
    gld = os.path.join(data_folder,'grilles_gld_xH/P1BP--2017030824934340CP_H1.hd')
    gri = loc.gld_xH(gld) 
    
   
    return mntbsq,gri
 
 
"""
calcul_gld3d -all -path_visee ./grilles_gld_xH -nom_visee P1BP--2017030824934340CP_H1 -type_visee LocalisateurGrille_Directe \
-mnt ./MNT_extrait/mnt_extrait -repter GRS80:G-D/:H-M -path_grille . -nom_grille test_intersect_euclide -convention BABEL \
-format BSQ -nbcol 200 -nblig 200 -pascol 50 -paslig 60 -j0 0 -i0 0 -col0 100.5 -lig0 100.5 -matrice 1 0 0 1
"""
 
@pytest.mark.unit_tests
def test_gld_mnt():  
    """
    Test loc direct grid on dtm function
    """
    mntbsq,gri = prepare_loc()                  
    gri_gld = gri.fct_gld_mnt(100.5, 100.5, 60, 50, 200, 200,mntbsq)                
    assert(True)                


@pytest.mark.unit_tests
def test_loc_dir_check_cube_mnt():  
    """
    Test direct localization check mnt cube
    """            
    assert(True)

@pytest.mark.unit_tests
def test_loc_dir_interp_visee_unitaire_gld():  
    """
    Test los interpolation
    """            
    assert(True)

    

@pytest.mark.unit_tests
def test_loc_dir_h():  
    """
    Test direct localization at constant altitude
    """
    lig = 100.5
    col = 50.5
    h = 100.0
    ___,gri = prepare_loc()
    lonlatalt = gri.fct_locdir_h(lig, col, h)
    # 57.2170045762374 21.959087150369 100
    valid_lon = 57.2170054518422
    valid_lat = 21.9590529453258
    valid_alt = 100.0
    diff_lon = lonlatalt[0] - valid_lon
    diff_lat = lonlatalt[1] - valid_lat
    diff_alt = lonlatalt[2] - valid_alt
    print("direct localization at constant altitude lig : {} col {} alt {}".format(lig,col,h))
    print("lon {} lat {} alt {} ".format(lonlatalt[0],lonlatalt[1],lonlatalt[2]))
    print('diff_lon {} diff_lat {} diff_alt {}'.format(diff_lon, diff_lat, diff_alt))
    assert(valid_lon == pytest.approx(lonlatalt[0],abs=1e-12))
    assert(valid_lat == pytest.approx(lonlatalt[1],abs=1e-12))
    assert(valid_alt == pytest.approx(lonlatalt[2],abs=1e-8))

@pytest.mark.unit_tests
def test_loc_dir_mnt():  
    """
    Test direct localization on DTM
    """            
    assert(True)

@pytest.mark.unit_tests
def test_loc_dir_mnt_opt():  
    """
    Test direct localization on DTM
    """            
    assert(True)

@pytest.mark.unit_tests
def test_loc_inv():
    """
    Test inverse localization
    """
    #init des predicteurs
    ___,gri = prepare_loc()
    gri.init_pred_loc_inv()
    lon = 57.2167252772905
    lat = 21.9587514585812
    alt = 10.0
    #init des predicteurs

    inv_lig,inv_col,valid = gri.fct_locinv([lon,lat,alt])
    valid_lig = 50.5
    valid_col = 10.0
    print("inverse localization  : lon {} lat {} alt {}".format(lon,lat,alt))
    print("lig {} col {}  ".format(inv_lig, inv_col))
    print('diff_lig {} diff_col {} '.format(inv_lig - valid_lig, inv_col - valid_col))   
    assert(inv_lig == pytest.approx(valid_lig,abs=1e-2))
    assert(inv_col == pytest.approx(valid_col,abs=1e-2))     

@pytest.mark.unit_tests
def test_loc_intersection():
    """
    Test direct localization intersection function
    """
    assert(True)



@pytest.mark.unit_tests
def test_loc_dir_loc_inv():
    """
    Test direct localization followed by inverse one
    """
    lig = 150.5
    col = 20.5
    h = 10.0
    ___,gri = prepare_loc()
    #init des predicteurs
    gri.init_pred_loc_inv()
    lonlatalt = gri.fct_locdir_h(lig, col, h)
    inv_lig,inv_col,valid = gri.fct_locinv(lonlatalt)

    print('lig {} col {} valid {}'.format(inv_lig, inv_col, valid))
    assert(lig == pytest.approx(inv_lig,abs=1e-2))
    assert(col == pytest.approx(inv_col,abs=1e-2))
    assert(valid == 1)


                
@pytest.mark.unit_tests
def test_coloc():
    """
    Test coloc function
    """
    mntbsq,gri = prepare_loc()
    gri.init_pred_loc_inv()
    gricol = loc.fct_coloc(gri, gri, mntbsq, 0.5, 0.5, 10, 100, 20, 20)
    assert(True)
    



    
    



