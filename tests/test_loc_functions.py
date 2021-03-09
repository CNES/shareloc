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

"""
Test module for localisation class shareloc/localisation.py
"""


import os
import pytest
import numpy as np
from utils import test_path

from shareloc.grid import Grid, coloc
from shareloc.dtm import DTM
from shareloc.rpc.rpc import RPC
from shareloc.localization import Localization
from shareloc.localization import coloc as coloc_rpc


def prepare_loc(alti="geoide", id_scene="P1BP--2017030824934340CP"):
    """
    Read multiH grid
    :param alti: alti validation dir
    :param id_scene: scene euclidium ID
    :return: multi H grid
    :rtype: str
    """
    data_folder = test_path(alti, id_scene)
    # chargement du mnt
    fic = os.path.join(data_folder, "MNT_extrait/mnt_extrait.c1")
    dtmbsq = DTM(fic)
    # chargement des grilles
    gld = os.path.join(data_folder, "grilles_gld_xH/{}_H1.hd".format(id_scene))
    gri = Grid(gld)
    return dtmbsq, gri


@pytest.mark.parametrize("idxcol,idxrow", [(20, 10)])
@pytest.mark.parametrize("row0,col0,steprow,stepcol,nbrow,nbcol", [(100.5, 150.5, 60, 50, 200, 200)])
@pytest.mark.unit_tests
def test_gld_dtm(idxrow, idxcol, row0, col0, steprow, stepcol, nbrow, nbcol):
    """
    Test loc direct grid on dtm function
    """
    dtmbsq, gri = prepare_loc()

    gri_gld = gri.direct_loc_grid_dtm(row0, col0, steprow, stepcol, nbrow, nbcol, dtmbsq)

    lonlatalt = gri_gld[:, idxrow, idxcol]

    row = row0 + steprow * idxrow
    col = col0 + stepcol * idxcol

    valid_lonlatalt = gri.direct_loc_dtm(row, col, dtmbsq)
    print("lon {} lat {} alt {} ".format(lonlatalt[0], lonlatalt[1], lonlatalt[2]))
    print("valid lon {} lat {} alt {} ".format(valid_lonlatalt[0], valid_lonlatalt[1], valid_lonlatalt[2]))
    assert lonlatalt == pytest.approx(valid_lonlatalt, abs=1e-12)


@pytest.mark.parametrize("col,row", [(50.5, 100.5)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.21700176041541, 21.959197148974, 238.0)])
@pytest.mark.unit_tests
def test_loc_dir_instersect_cube_dtm(col, row, valid_lon, valid_lat, valid_alt):
    """
    Test direct localization check dtm cube
    """
    dtmbsq, gri = prepare_loc()

    visee = np.zeros((3, gri.nbalt))
    vislonlat = gri.interpolate_grid_in_plani(row, col)
    visee[0, :] = vislonlat[0]
    visee[1, :] = vislonlat[1]
    visee[2, :] = gri.alts_down
    visee = visee.T
    (__, __, position, __) = dtmbsq.intersect_dtm_cube(visee)
    assert position[0] == pytest.approx(valid_lon, abs=1e-12)
    assert position[1] == pytest.approx(valid_lat, abs=1e-12)
    assert position[2] == pytest.approx(valid_alt, abs=1e-12)


@pytest.mark.parametrize("col,row", [(50.5, 100.5)])
@pytest.mark.parametrize("valid_lon,valid_lat", [(57.2169100597702, 21.96277930762832)])
@pytest.mark.unit_tests
def test_loc_dir_interp_visee_unitaire_gld(row, col, valid_lon, valid_lat):
    """
    Test los interpolation
    """
    ___, gri = prepare_loc()
    visee = gri.interpolate_grid_in_plani(row, col)
    print(visee)
    assert visee[0][1] == pytest.approx(valid_lon, abs=1e-12)
    assert visee[1][1] == pytest.approx(valid_lat, abs=1e-12)


@pytest.mark.parametrize("col,row,h", [(50.5, 100.5, 100.0)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.2170054518422, 21.9590529453258, 100.0)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_h(col, row, h, valid_lon, valid_lat, valid_alt):
    """
    Test direct localization at constant altitude
    """
    ___, gri = prepare_loc()
    loc = Localization(gri, elevation=None)
    lonlatalt = loc.direct(row, col, h)

    diff_lon = lonlatalt[0] - valid_lon
    diff_lat = lonlatalt[1] - valid_lat
    diff_alt = lonlatalt[2] - valid_alt
    print("direct localization at constant altitude row : {} col {} alt {}".format(row, col, h))
    print("lon {} lat {} alt {} ".format(lonlatalt[0], lonlatalt[1], lonlatalt[2]))
    print("diff_lon {} diff_lat {} diff_alt {}".format(diff_lon, diff_lat, diff_alt))
    assert valid_lon == pytest.approx(lonlatalt[0], abs=1e-12)
    assert valid_lat == pytest.approx(lonlatalt[1], abs=1e-12)
    assert valid_alt == pytest.approx(lonlatalt[2], abs=1e-8)


@pytest.mark.parametrize("col,row,h", [(150.5, 100, 100.0)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_vs_loc_rpc(row, col, h):
    """
    Test direct localization coherence with grid and RPC
    """
    id_scene = "P1BP--2018122638935449CP"
    ___, gri = prepare_loc("ellipsoide", id_scene)
    loc_grid = Localization(gri)
    # init des predicteurs
    lonlatalt = loc_grid.direct(row, col, h)

    data_folder = test_path()
    fichier_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_any(fichier_dimap)

    loc_rpc = Localization(fctrat)
    # grid (from physical model) and RP have 0.5 pixel shift
    lonlatalt_rpc = loc_rpc.direct(row - 0.5, col - 0.5, h)

    diff_lon = lonlatalt[0] - lonlatalt_rpc[0]
    diff_lat = lonlatalt[1] - lonlatalt_rpc[1]
    diff_alt = lonlatalt[2] - lonlatalt_rpc[2]
    assert diff_lon == pytest.approx(0.0, abs=1e-7)
    assert diff_lat == pytest.approx(0.0, abs=1e-7)
    assert diff_alt == pytest.approx(0.0, abs=1e-7)


@pytest.mark.parametrize("index_x,index_y", [(10.5, 20.5)])
@pytest.mark.unit_tests
def test_sensor_loc_dir_dtm(index_x, index_y):
    """
    Test direct localization on DTM
    """
    dtmbsq, gri = prepare_loc()
    loc = Localization(gri, elevation=dtmbsq)

    vect_index = [index_x, index_y]
    [lon, lat] = dtmbsq.index_to_ter(vect_index)
    print([lon, lat])
    alt = dtmbsq.interpolate(index_x - 0.5, index_y - 0.5)

    row, col, alt = loc.inverse(lon, lat, alt)
    print("row col ", row, col)
    lonlath = loc.direct(row, col)
    assert lon == pytest.approx(lonlath[0], abs=1e-8)
    assert lat == pytest.approx(lonlath[1], abs=1e-8)
    assert alt == pytest.approx(lonlath[2], abs=1e-4)


@pytest.mark.parametrize("valid_row,valid_col", [(50.5, 10.0)])
@pytest.mark.parametrize("lon,lat,alt", [(57.2167252772905, 21.9587514585812, 10.0)])
@pytest.mark.unit_tests
def test_sensor_loc_inv(lon, lat, alt, valid_col, valid_row):
    """
    Test inverse localization
    """

    ___, gri = prepare_loc()

    loc = Localization(gri)
    inv_row, inv_col, h = loc.inverse(lon, lat, alt)

    print("inverse localization  : lon {} lat {} alt {}".format(lon, lat, alt))
    print("row {} col {}  ".format(inv_row, inv_col))
    print("diff_row {} diff_col {} ".format(inv_row - valid_row, inv_col - valid_col))
    assert inv_row == pytest.approx(valid_row, abs=1e-2)
    assert inv_col == pytest.approx(valid_col, abs=1e-2)
    assert h == alt


@pytest.mark.parametrize("lon,lat,alt", [(2.12026631, 31.11245154, 10.0)])
@pytest.mark.unit_tests
def test_sensor_loc_inv_vs_loc_rpc(lon, lat, alt):
    """
    Test direct localization coherence with grid and RPC
    """
    id_scene = "P1BP--2018122638935449CP"
    ___, gri = prepare_loc("ellipsoide", id_scene)
    loc_grid = Localization(gri)

    [row, col, __] = loc_grid.inverse(lon, lat, alt)
    data_folder = test_path()
    fichier_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))

    fctrat = RPC.from_any(fichier_dimap, topleftconvention=True)

    loc_rpc = Localization(fctrat)
    [row_rpc, col_rpc, __] = loc_rpc.inverse(lon, lat, alt)
    diff_row = row_rpc - row
    diff_col = col_rpc - col
    # delta vt 0.5 pixel shift between physical model and rpc OTB
    assert diff_row == pytest.approx(-0.5, abs=1e-2)
    assert diff_col == pytest.approx(-0.5, abs=1e-2)


@pytest.mark.parametrize(
    "col_min_valid",
    [([4.15161251e-02, 1.95057636e-01, 1.10977819e00, -8.35016563e-04, -3.50772271e-02, -9.46432481e-03])],
)
@pytest.mark.parametrize(
    "row_min_valid", [([0.05440845, 1.26513831, -0.36737151, -0.00229532, -0.07459378, -0.02558954])]
)
@pytest.mark.parametrize(
    "col_max_valid",
    [([1.76451389e-02, 2.05533045e-01, 1.11758291e00, -9.50086076e-04, -3.59923603e-02, -1.03291594e-02])],
)
@pytest.mark.parametrize(
    "row_max_valid", [([0.07565692, 1.27499912, -0.36677813, -0.00252395, -0.07539624, -0.0270914])]
)
@pytest.mark.parametrize("valid_offset_lon", [([57.37295223744326, 0.15660032225072484])])
@pytest.mark.parametrize("valid_offset_lat", [([22.066877016445275, 0.14641205050748773])])
@pytest.mark.parametrize("valid_offset_row", [([24913.0, 24912.5])])
@pytest.mark.parametrize("valid_offset_col", [([19975.5, 19975.0])])
@pytest.mark.unit_tests
def test_pred_loc_inv(
    col_min_valid,
    row_min_valid,
    col_max_valid,
    row_max_valid,
    valid_offset_lon,
    valid_offset_lat,
    valid_offset_row,
    valid_offset_col,
):
    """
    Test inverse localization
    """
    # init des predicteurs
    ___, gri = prepare_loc()
    gri.estimate_inverse_loc_predictor()

    assert gri.pred_col_min.flatten() == pytest.approx(col_min_valid, abs=1e-6)
    assert gri.pred_row_min.flatten() == pytest.approx(row_min_valid, abs=1e-6)
    assert gri.pred_col_max.flatten() == pytest.approx(col_max_valid, abs=1e-6)
    assert gri.pred_row_max.flatten() == pytest.approx(row_max_valid, abs=1e-6)
    assert gri.pred_ofset_scale_lon == pytest.approx(valid_offset_lon, abs=1e-12)
    assert gri.pred_ofset_scale_lat == pytest.approx(valid_offset_lat, abs=1e-12)
    assert gri.pred_ofset_scale_row == pytest.approx(valid_offset_row, abs=1e-6)
    assert gri.pred_ofset_scale_col == pytest.approx(valid_offset_col, abs=1e-6)


@pytest.mark.parametrize("col,row", [(50.5, 100.5)])
@pytest.mark.parametrize("valid_lon,valid_lat,valid_alt", [(57.21700367698209, 21.95912227930429, 166.351227229112)])
@pytest.mark.unit_tests
def test_loc_intersection(row, col, valid_lon, valid_lat, valid_alt):
    """
    Test direct localization intersection function
    """
    dtmbsq, gri = prepare_loc()

    visee = np.zeros((3, gri.nbalt))
    vislonlat = gri.interpolate_grid_in_plani(row, col)
    visee[0, :] = vislonlat[0]
    visee[1, :] = vislonlat[1]
    visee[2, :] = gri.alts_down
    visee = visee.T
    (__, __, point_b, alti) = dtmbsq.intersect_dtm_cube(visee)
    (__, __, point_dtm) = dtmbsq.intersection(visee, point_b, alti)
    assert point_dtm[0] == pytest.approx(valid_lon, abs=1e-12)
    assert point_dtm[1] == pytest.approx(valid_lat, abs=1e-12)
    assert point_dtm[2] == pytest.approx(valid_alt, abs=1e-12)


@pytest.mark.parametrize("col,row,h", [(20.5, 150.5, 10.0)])
@pytest.mark.unit_tests
def test_loc_dir_loc_inv(row, col, h):
    """
    Test direct localization followed by inverse one
    """
    ___, gri = prepare_loc()
    # init des predicteurs
    gri.estimate_inverse_loc_predictor()
    (lon, lat, alt) = gri.direct_loc_h(row, col, h)
    inv_row, inv_col, h = gri.inverse_loc(lon, lat, alt)

    print("row {} col {} ".format(inv_row, inv_col))
    assert row == pytest.approx(inv_row, abs=1e-2)
    assert col == pytest.approx(inv_col, abs=1e-2)


# delta vt 0.5 pixel shift between physical model and rpc OTB
@pytest.mark.parametrize(
    "id_scene, rpc, col,row, h, delta_vt",
    [
        ("P1BP--2018122638935449CP", "PHRDIMAP_P1BP--2018122638935449CP.XML", 150.5, 20.5, 10.0, 0.5),
        ("P1BP--2017092838284574CP", "RPC_PHR1B_P_201709281038045_SEN_PRG_FC_178608-001.XML", 150.5, 20.5, 10.0, 0.5),
    ],
)
@pytest.mark.unit_tests
def test_loc_dir_loc_inv_rpc(id_scene, rpc, row, col, h, delta_vt):
    """
    Test direct localization followed by inverse one
    """
    ___, gri = prepare_loc("ellipsoide", id_scene)
    # init des predicteurs
    gri.estimate_inverse_loc_predictor()
    lonlatalt = gri.direct_loc_h(row, col, h)

    data_folder = test_path()
    fichier_dimap = os.path.join(data_folder, "rpc", rpc)

    fctrat = RPC.from_any(fichier_dimap, topleftconvention=True)
    (inv_row, inv_col, __) = fctrat.inverse_loc(lonlatalt[0], lonlatalt[1], lonlatalt[2])
    print("row {} col {}".format(inv_row, inv_col))

    assert row == pytest.approx(inv_row + delta_vt, abs=1e-2)
    assert col == pytest.approx(inv_col + delta_vt, abs=1e-2)


@pytest.mark.parametrize("l0_src,c0_src, steprow_src, stepcol_src,nbrow_src,nbcol_src", [(0.5, 1.5, 10, 100, 20, 20)])
@pytest.mark.parametrize("col,row", [(1, 3)])
@pytest.mark.unit_tests
def test_coloc(l0_src, c0_src, steprow_src, stepcol_src, nbrow_src, nbcol_src, row, col):
    """
    Test coloc function
    """
    dtmbsq, gri = prepare_loc()
    gri.estimate_inverse_loc_predictor()

    gricol = coloc(gri, gri, dtmbsq, [l0_src, c0_src], [steprow_src, stepcol_src], [nbrow_src, nbcol_src])

    assert gricol[0, row, col] == pytest.approx(row * steprow_src + l0_src, 1e-6)
    assert gricol[1, row, col] == pytest.approx(col * stepcol_src + c0_src, 1e-6)


@pytest.mark.parametrize("index_x,index_y", [(10.5, 20.5)])
@pytest.mark.parametrize("valid_alt", [(198.0)])
@pytest.mark.unit_tests
def test_interp_dtm(index_x, index_y, valid_alt):
    """
    Test coloc function
    """
    dtmbsq, ___ = prepare_loc()

    vect_index = [index_x, index_y]
    coords = dtmbsq.index_to_ter(vect_index)
    print(coords)
    alti = dtmbsq.interpolate(index_x - 0.5, index_y - 0.5)
    assert alti == valid_alt


@pytest.mark.parametrize("col,lig,h", [(1000.5, 1500.5, 10.0)])
@pytest.mark.unit_tests
def test_loc_dir_loc_inv_couple(lig, col, h):
    """
    Test direct localization followed by inverse one for phr couple
    """
    id_scene_right = "P1BP--2017092838319324CP"
    ___, gri_right = prepare_loc("ellipsoide", id_scene_right)
    id_scene_left = "P1BP--2017092838284574CP"
    ___, gri_left = prepare_loc("ellipsoide", id_scene_left)
    # init des predicteurs
    gri_right.estimate_inverse_loc_predictor()
    lonlatalt = gri_left.direct_loc_h(lig, col, h)
    inv_lig, inv_col, __ = gri_right.inverse_loc(lonlatalt[0], lonlatalt[1], lonlatalt[2])

    print("lig {} col {}".format(inv_lig, inv_col))
    # assert(lig == pytest.approx(inv_lig,abs=1e-2))
    # assert(col == pytest.approx(inv_col,abs=1e-2))
    # assert(valid == 1)


@pytest.mark.parametrize("col,row,alt", [(600, 200, 125)])
def test_colocalization(col, row, alt):
    """
    Test colocalization function using rpc
    """

    data_folder = test_path()
    id_scene = "P1BP--2018122638935449CP"
    file_dimap = os.path.join(data_folder, "rpc/PHRDIMAP_{}.XML".format(id_scene))
    fctrat = RPC.from_dimap_v1(file_dimap)

    row_coloc, col_coloc, _ = coloc_rpc(fctrat, fctrat, row, col, alt)

    assert row == pytest.approx(row_coloc, abs=1e-1)
    assert col == pytest.approx(col_coloc, abs=1e-1)
