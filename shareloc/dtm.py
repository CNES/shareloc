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
This module contains the DTM class to handle dtm intersection.
DTM files must be in BSQ format.
"""

import numpy as np
from shareloc.image.dtm_image import DTMImage
from shareloc.math_utils import interpol_bilin


class DTM:
    """DTM class dedicated
    works only with BSQ file format
    we work in cell convention [0,0] is the first cell center (not [0.5,0.5])
    """

    # gitlab issue #56
    # pylint: disable=too-many-instance-attributes
    def __init__(self, dtm_filename, dtm_format="bsq"):
        """
        Constructor
        :param dtm_filename: dtm filename
        :type dtm_filename: string
        :param dtm_format: grid format (by default bsq)
        :type dtm_format: string
        """
        self.dtm_file = dtm_filename
        self.format = dtm_format
        self.alt_data = None
        self.alt_min = None
        self.alt_max = None
        self.origin_x = None
        self.origin_y = None
        self.pixel_size_x = None
        self.pixel_size_y = None
        self.column_nb = None
        self.row_nb = None
        self.plane_coef_a = None
        self.plane_coef_b = None
        self.plane_coef_c = None
        self.plane_coef_d = None
        self.alt_min_cell = None
        self.alt_max_cell = None
        self.tol_z = 0.0001

        # lecture mnt
        self.dtm_image = DTMImage(self.dtm_file, read_data=True)
        self.alt_data = self.dtm_image.data

        self.init_min_max()
        self.alt_max = self.alt_data.max()
        self.alt_min = self.alt_data.min()

        self.plane_coef_a = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.plane_coef_b = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        self.plane_coef_c = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        self.plane_coef_d = np.array(
            [0.0, self.dtm_image.nb_rows - 1.0, 0.0, self.dtm_image.nb_cols - 1.0, self.alt_min, self.alt_max]
        )

        self.plans = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, self.dtm_image.nb_rows - 1.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, self.dtm_image.nb_cols - 1.0],
                [0.0, 0.0, 1.0, self.alt_min],
                [0.0, 0.0, 1.0, self.alt_max],
            ]
        )

    def eq_plan(self, i, position):
        """
        return evaluation of equation on a plane on DTM cube
        :param i: face index
        :type i: int
        :param position: position
        :type position: numpy.array (1x3)
        :return evaluation on the plan
        :rtype float
        """
        return (
            self.plane_coef_a[i] * position[0]
            + self.plane_coef_b[i] * position[1]
            + self.plane_coef_c[i] * position[2]
            - self.plane_coef_d[i]
        )

    def ter_to_index(self, vect_ter):
        """
        terrain to index conversion
        :param vect_ter: terrain coordinate
        :type vect_ter: array (1x2 or 1x3) if dimension is 3 , last coordinate is unchanged (alt)
        :return index coordinates
        :rtype array (1x2 or 1x3)
        """
        vect_dtm = vect_ter.copy()
        (vect_dtm[0], vect_dtm[1]) = self.dtm_image.transform_physical_point_to_index(vect_ter[1], vect_ter[0])
        return vect_dtm

    def ters_to_indexs(self, vect_ters):
        """
        terrain to index conversion
        :param vect_ters: terrain coordinates
        :type vect_ters: array (nx2 or nx3) if dimension is 3 , last coordinates is iuchanged (alt)
        :return index coordinates
        :rtype array (nx2 or nx3)
        """
        vect_dtms = vect_ters.copy()
        for i, vect_ter in enumerate(vect_ters):
            vect_dtms[i, :] = self.ter_to_index(vect_ter)
        return vect_dtms

    def index_to_ter(self, vect_dtm):
        """
        index to terrain conversion
        :param vect_dtm: index coordinate
        :type vect_dtm: array (1x2 or 1x3) if dimension is 3 , last coordinate is unchanged (alt)
        :return terrain coordinates
        :rtype array (1x2 or 1x3)
        """
        vect_ter = vect_dtm.copy()
        (vect_ter[1], vect_ter[0]) = self.dtm_image.transform_index_to_physical_point(vect_dtm[0], vect_dtm[1])
        return vect_ter

    def interpolate(self, pos_row, pos_col):
        """
        interpolate altitude
        if interpolation is done outside DTM the penultimate index is used (or first if d is negative).
        :param pos_row: cell position row
        :type pos_row: float
        :param pos_col: cell position col
        :type pos_col: float
        :return interpolated altitude
        :rtype float
        """
        alt = interpol_bilin(
            [self.alt_data[np.newaxis, :, :]], self.dtm_image.nb_rows, self.dtm_image.nb_cols, pos_row, pos_col
        )[0][0]
        return alt

    def init_min_max(self):
        """
        initialize min/max at each dtm cell
        """
        column_nb_minus = self.dtm_image.nb_cols - 1
        row_nb_minus = self.dtm_image.nb_rows - 1

        self.alt_max_cell = np.zeros((row_nb_minus, column_nb_minus))
        self.alt_min_cell = np.zeros((row_nb_minus, column_nb_minus))

        # On calcule les altitudes min et max du MNT
        n_min = 32000
        n_max = -32000

        for i in range(row_nb_minus):
            k = 0
            for j in range(column_nb_minus):
                k += 1
                d_alt_min = n_min
                d_alt_max = n_max

                d_z1 = self.alt_data[i, j]
                if d_z1 < d_alt_min:
                    d_alt_min = d_z1
                if d_z1 > d_alt_max:
                    d_alt_max = d_z1

                d_z2 = self.alt_data[i, j + 1]
                if d_z2 < d_alt_min:
                    d_alt_min = d_z2
                if d_z2 > d_alt_max:
                    d_alt_max = d_z2

                d_z3 = self.alt_data[i + 1, j]
                if d_z3 < d_alt_min:
                    d_alt_min = d_z3
                if d_z3 > d_alt_max:
                    d_alt_max = d_z3

                d_z4 = self.alt_data[i + 1, j + 1]
                if d_z4 < d_alt_min:
                    d_alt_min = d_z4
                if d_z4 > d_alt_max:
                    d_alt_max = d_z4

                # GDN Correction BUG Interesctor
                # Il ne faut surtout pas prendre l'arrondi pour plusieurs raisons
                # 1. l'algo ulterieur ne resiste pas touours bien lorsque la maille est plate,
                #    il est donc deconseille de fournir en sortie i_altmin = i_altmax
                #    a moins que ce soi vraiment le cas en valeurs reelles
                # 2. il ne faut pas initialiser les mailles consecutives de la meme maniere par arrondi
                #    car si l'altitude min de l'une correspond a l'altitude max de l'autre il faut les distinguer
                #    par un ceil et un floor pour que les cubes se chevauchent legerement en altitude
                #    et non pas jointifs strictement
                i_altmin = np.floor(d_alt_min)
                i_altmax = np.ceil(d_alt_max)
                self.alt_min_cell[i, j] = i_altmin
                self.alt_max_cell[i, j] = i_altmax

    # gitlab issue #56
    # pylint: disable=too-many-branches
    def intersect_dtm_cube(self, los):
        """
        DTM cube intersection
        :param los :  line of sight
        :type los : numpy.array
        :return intersection information (True,an intersection has been found ?, (lon,lat) of dtm position, altitude)
        :rtype tuple (bool, bool, numpy.array, float)
        """
        # los: (n,3):
        point_b = None
        h_intersect = None
        (coord_col_i, coord_row_i, coord_alt_i, alti_layer_i) = ([], [], [], [])
        nbalt = los.shape[0]
        los_index = self.ters_to_indexs(los)
        # -----------------------------------------------------------------------
        # Nombre d'intersections valides trouvees
        nbi = 0
        # -----------------------------------------------------------------------
        # On boucle sur les plans du cube DTM
        for plane_index in range(6):
            # -----------------------------------------------------------------------
            # Initialisation du sommet de la visee
            los_hat = los_index[0, :]
            # -----------------------------------------------------------------------
            # Initialisation de la position par / au plan
            # print self.plans[plane_index,:-1]
            # los_hat_onplane = (self.plans[plane_index,:-1]*los_hat).sum() - self.d[plane_index]
            # print los_hat_onplane
            los_hat_onplane = self.eq_plan(plane_index, los_hat)
            # -----------------------------------------------------------------------
            # On boucle sur les segments de la visee et on controle
            # si on traverse ou pas la face courante plane_index du cube DTM
            for alti_layer in range(nbalt):
                # -----------------------------------------------------------------------
                # Transfert du point B dans le point A
                los_a = los_hat_onplane
                s_a = los_hat.copy()
                # -----------------------------------------------------------------------
                # Reinit du point B
                los_hat = los_index[alti_layer, :]  # dg on itere sur les differents points de la visee
                # print los_hat
                # -----------------------------------------------------------------------
                # Initialisation de la position par / au plan
                los_hat_onplane = self.eq_plan(plane_index, los_hat)
                # print 'posAposB',los_a,los_hat_onplane
                # -----------------------------------------------------------------------
                # Test d'intersection : los_a et los_hat_onplane de signes opposes
                if los_a * los_hat_onplane <= 0:
                    if not los_a and not los_hat_onplane:
                        # Trop de solutions !! (A et B sont sur le plan) #comment on le gere ??
                        print("too many solutions")
                    # -----------------------------------------------------------------------
                    if not los_a:
                        # A est solution (il est sur le plan)
                        coord_col_i.append(s_a[0])
                        coord_row_i.append(s_a[1])
                        coord_alt_i.append(s_a[2])
                        alti_layer_i.append(alti_layer - 1)
                    # -----------------------------------------------------------------------
                    elif not los_hat_onplane:
                        # B est solution (il est sur le plan)
                        coord_col_i.append(los_hat[0])
                        coord_row_i.append(los_hat[1])
                        coord_alt_i.append(los_hat[2])
                        alti_layer_i.append(alti_layer)
                    # -----------------------------------------------------------------------
                    else:
                        # -----------------------------------------------------------------------
                        # A et B sont de part et d'autre du plan
                        # Coefficients d'interpolation de l'intersection
                        # entre A et B
                        interp_coef_a = los_hat_onplane / (los_hat_onplane - los_a)
                        interp_coef_b = -los_a / (los_hat_onplane - los_a)
                        # Affectation ou interpolation
                        # NB : pour eviter les pb lors du test
                        #      <estSurCube> (voir + loin)
                        # . coordonn?e <u> (ligne)
                        # -----------------------------------------------------------------------
                        if plane_index < 2:
                            coord_col_i.append(self.plane_coef_d[plane_index])
                        # -----------------------------------------------------------------------
                        else:
                            coord_col_i.append(interp_coef_a * s_a[0] + interp_coef_b * los_hat[0])
                        # -----------------------------------------------------------------------
                        # . coordonnee <v> (colonne)
                        if 1 < plane_index < 4:
                            coord_row_i.append(self.plane_coef_d[plane_index])
                        else:
                            coord_row_i.append(interp_coef_a * s_a[1] + interp_coef_b * los_hat[1])
                        # -----------------------------------------------------------------------
                        # . coordonn?e <z> (altitude)
                        if plane_index > 3:
                            coord_alt_i.append(self.plane_coef_d[plane_index])
                        # -----------------------------------------------------------------------
                        else:
                            coord_alt_i.append(interp_coef_a * s_a[2] + interp_coef_b * los_hat[2])
                        # . coordonn?e <h> (abscisse visee)
                        alti_layer_i.append(alti_layer - interp_coef_a)  # index non entier de l'intersection
                    # -----------------------------------------------------------------------
                    # Incrementation du nombre d'intersections trouvees
                    nbi += 1
                    # -----------------------------------------------------------------------
                    # Passage a la face du cube suivante
                    break

        # -----------------------------------------------------------------------
        # Tri des points le long de la visee (il y en a au moins deux)
        # on les range par ordre decroissant en altitude
        for alti_layer in range(nbi):
            for next_alti_layer in range(alti_layer + 1, nbi):
                if alti_layer_i[next_alti_layer] < alti_layer_i[alti_layer]:
                    dtmp = coord_col_i[alti_layer]
                    coord_col_i[alti_layer] = coord_col_i[next_alti_layer]
                    coord_col_i[next_alti_layer] = dtmp
                    dtmp = coord_row_i[alti_layer]
                    coord_row_i[alti_layer] = coord_row_i[next_alti_layer]
                    coord_row_i[next_alti_layer] = dtmp
                    dtmp = coord_alt_i[alti_layer]
                    coord_alt_i[alti_layer] = coord_alt_i[next_alti_layer]
                    coord_alt_i[next_alti_layer] = dtmp
                    dtmp = alti_layer_i[alti_layer]
                    alti_layer_i[alti_layer] = alti_layer_i[next_alti_layer]
                    alti_layer_i[next_alti_layer] = dtmp

        # -----------------------------------------------------------------------
        # Filtrage des points non situes sur le cube
        alti_layer = 0
        while alti_layer < nbi:
            # test a l'interieur du cube
            test_on_cube = (
                (coord_col_i[alti_layer] >= self.plane_coef_d[0])
                and (coord_col_i[alti_layer] <= self.plane_coef_d[1])
                and (coord_row_i[alti_layer] >= self.plane_coef_d[2])
                and (coord_row_i[alti_layer] <= self.plane_coef_d[3])
                and (coord_alt_i[alti_layer] >= self.plane_coef_d[4])
                and (coord_alt_i[alti_layer] <= self.plane_coef_d[5])
            )
            if not test_on_cube:
                # On translate tous les points suivants (on ecrase ce point non valide)
                for next_alti_layer in range(alti_layer + 1, nbi):
                    coord_col_i[next_alti_layer - 1] = coord_col_i[next_alti_layer]
                    coord_row_i[next_alti_layer - 1] = coord_row_i[next_alti_layer]
                    coord_alt_i[next_alti_layer - 1] = coord_alt_i[next_alti_layer]
                    alti_layer_i[next_alti_layer - 1] = alti_layer_i[next_alti_layer]
                nbi -= 1
            else:
                alti_layer += 1
        # -----------------------------------------------------------------------
        # Pas de solution si 0 ou 1 seul point trouve (on a tangente le cube)
        if nbi < 2:
            b_trouve = False
            return True, b_trouve, point_b, h_intersect
        # -----------------------------------------------------------------------
        # Il ne reste que 2 points donc on traverse le cube
        # LAIG-FA-MAJA-2168-CNES: plus de filtrage sur les point identiques. Il peut y avoir un nombre depoints > 2
        # Initialisation du point courant
        # Coordonnees DTM
        point_dtm = np.zeros(3)
        point_dtm[0] = coord_col_i[0]
        point_dtm[1] = coord_row_i[0]
        point_dtm[2] = coord_alt_i[0]
        # PointDTM est la premiere intersection avec le cube (lig, col)
        # -----------------------------------------------------------------------
        # h dans gld 3D
        h_intersect = alti_layer_i[0]
        # h_intersect correspond a l'index (non entier) d'interpolation en h
        # -----------------------------------------------------------------------
        # Coordonnees terrain
        point_b = self.index_to_ter(point_dtm)
        # point_b est le point Terrain (lon,lat)
        # -----------------------------------------------------------------------
        # Fin, retour
        b_trouve = True
        return (True, b_trouve, point_b, h_intersect)

    # gitlab issue #56
    # pylint: disable=too-many-locals
    # pylint: disable=too-many-function-args
    # pylint: disable=too-many-nested-blocks
    # pylint: disable=too-many-statements
    def intersection(self, los, point_b, h_intersect):
        """
        DTM intersection
        :param los :  line of sight
        :type los : numpy.array
        :param point_b :  position of intersection in DTM cube
        :type point_b : numpy.array
        :param h_intersect :  altitude in DTM cube
        :type h_intersect : float
        :return intersection information (True,an intersection has been found ?, position of intersection)
        :rtype tuple (bool, bool, numpy.array)
        """
        los_index = self.ters_to_indexs(los)
        point_b_dtm = self.ter_to_index(point_b)
        point_r = np.zeros(3)
        (npl, _) = los.shape
        alti = range(npl, -1, -1)

        p_1 = point_b_dtm.copy()  # [p_1[0],p_1[1],p_1[2]]

        h_intersect_p1 = h_intersect

        n_row = self.dtm_image.nb_rows
        n_col = self.dtm_image.nb_cols

        # 1 - Initilialisation et tests prealables
        #   1.1 - Test si le sommet est au-dessu du DTM
        #       - Calcul de l'altitude du DTM ? la position du sommet
        alti_1 = self.interpolate(p_1[0], p_1[1])
        #       - Calcul de l'ecart d'altitude au DTM
        d_alti_1 = p_1[2] - alti_1

        #       - Test si le nouveau point haut est au dessus du DTM
        if d_alti_1 < 0:
            #       - Point situ? en dessous du DTM
            #          . ceci signifie que la vis?e rentre dans le DTM par le c?t?
            #          . donc en dessous, pas de solution
            b_trouve = False
            return True, b_trouve, point_r

        #   1.2 - Initialisation du rang du premier sommet de la visee
        i_0 = int(np.floor(h_intersect_p1))

        #   1.3 - Initialisation du point de depart (dans p_2)
        p_2 = point_b_dtm.copy()
        h_intersect_p2 = h_intersect

        # 2. - Boucle sur les plans de grille
        while i_0 < (los_index.size - 1):
            # 2.1 - Initialisation du sommet courant de la visee
            col_0 = los_index[i_0][0]
            row_0 = los_index[i_0][1]
            z_0 = los_index[i_0][2]
            z_1 = los_index[i_0 + 1][2]

            # 2.2 - Initialisation de la visee DTM
            los_dtm = los_index[i_0 + 1] - los_index[i_0]

            # 2.3 - Test si visee verticale
            if los_dtm[0] == 0 and los_dtm[1] == 0:
                # 2.3.1 - La vis?e est verticale :
                #    - Calcul de l'altitude du DTM ? la position du sommet
                alti_1 = self.interpolate(col_0, row_0)

                #    Test si le plan suivant est en dessous du DTM
                if los_index[i_0 + 1][2] <= alti_1:
                    # Init point de sortie
                    p_1[0] = col_0
                    p_1[1] = row_0
                    p_1[2] = alti_1
                    b_trouve = True
                    self.index_to_ter(p_1, point_r)
                    return True, b_trouve, point_r
                # Positionnement sur le sommet suivant
                i_0 += 1
            else:
                # 2.3.2 - La visee n'est pas verticale :
                #         elle v_a donc survoler le DTM
                #         . on peut donc poursuivre
                #         . reste ? d?montrer que la vis?e se crashe sur le DTM...
                #
                # Initialisation du point de d?part
                # Initialisation de son abscisse sur la vis?e
                a_2 = h_intersect_p2 - i_0

                # Correction d'un bug FA DG 10, a la reinitialisation a_2 peut valoir 1
                # Ce qui a pour consequence de sauter au segment suivant de la visee
                # Ainsi, on peut perdre des points sur une tranche d'altitude
                # Pour etre sur de scanner le segment de visee, on reinitialise
                # le point au depart du segment, soit a_2 = 0
                if a_2 >= 1.0:
                    a_2 = 0.0

                # Initialisation de la premi?re maille DTM intersect?e
                #  - Initialisation des indices de la maille
                col_c = int(np.floor(p_2[0]))
                row_c = int(np.floor(p_2[1]))

                # NB :    pr?caution avant de d?marrer :
                #        . on se met du bon cote de la maille
                #        . en principe, on ne doit pas sortir du DTM
                # On rentre par le bas, la maille DTM est la precedente
                if (p_2[0] == col_c) and (los_dtm[0] < 0):
                    col_c -= 1

                # On rentre par la gauche, la maille DTM est la precedente
                if (p_2[1] == row_c) and (los_dtm[1] < 0):
                    row_c -= 1

                # LDD - On est deja en dehors des limites, on s'arrete
                if not ((a_2 < 1) and -1 < col_c < (n_row - 1) and -1 < row_c < (n_col - 1)):
                    b_trouve = False
                    return True, b_trouve, point_r

                # Boucle de recherche iterative de la maille intersectee
                while (a_2 < 1) and -1 < col_c < (n_row - 1) and -1 < row_c < (n_col - 1):
                    # - Altitudes min et max de la maille
                    h_i = self.alt_min_cell[col_c, row_c]
                    h_s = self.alt_max_cell[col_c, row_c]

                    # - Transfert : le point bas devient le point haut
                    # a1 = a_2;
                    # p_1 devient p_2
                    p_1 = p_2.copy()
                    h_intersect_p1 = h_intersect_p2.copy()

                    # 4.2 - D?termination d'un nouveau point bas
                    #      - Test d'orientation de la vis?e
                    if not los_dtm[0]:
                        # 4.2.1 - La vis?e est orientee pile poil est-ouest
                        #   p_2[0] = p_1[0] ; // inutile, est deja initialise
                        #       - Test d'orientation de la visee
                        if los_dtm[1] < 0:
                            # 4.2.1.1 - La vis?e part plein ouest
                            p_2[1] = row_c
                            row_c -= 1
                        else:
                            # 4.2.1.2 - La vis?e part plein est
                            row_c += 1
                            p_2[1] = row_c

                        a_2 = (p_2[1] - row_0) / los_dtm[1]
                        p_2[2] = z_0 + a_2 * los_dtm[2]

                    elif not los_dtm[1]:
                        # 4.2.2 - La vis?e est orient?e nord-sud
                        #  p_2[1] = p_1[1] ;
                        #       - Test d'orientation de la visee
                        if los_dtm[0] < 0:
                            # 4.2.2.1 - La vis?e part plein nord
                            p_2[0] = col_c
                            col_c -= 1
                        else:
                            # 4.2.2.2 - La vis?e part plein sud
                            col_c += 1
                            p_2[0] = col_c

                        a_2 = (p_2[0] - col_0) / los_dtm[0]
                        p_2[2] = z_0 + a_2 * los_dtm[2]
                    else:
                        # 4.2.3 - La vis?e est quelconque
                        #            - D?termination du cot? de sortie
                        if (los_dtm[0] < 0) and (los_dtm[0] <= los_dtm[1]) and (los_dtm[0] <= -los_dtm[1]):
                            # 4.2.3.1 - Vis?e principalement orient?e nord
                            #             - Intersection avec le c?t? nord
                            a_2 = (col_c - col_0) / los_dtm[0]
                            p_2[1] = row_0 + a_2 * los_dtm[1]

                            if (p_2[1] > row_c) and (p_2[1] < (row_c + 1)):
                                # La vis?e sort par le nord
                                p_2[0] = col_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[1] < row_c:
                                # La vis?e sort par l'ouest
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[1] > (row_c + 1):
                                # La vis?e sort par l'est
                                row_c += 1
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[1] == row_c:
                                # La vis?e sort par le coin nord-ouest
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1
                                row_c -= 1

                            elif p_2[1] == (row_c + 1):
                                # La vis?e sort par le coin nord-est
                                p_2[0] = col_c
                                row_c += 1
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                        elif (los_dtm[1] > 0) and (los_dtm[1] >= los_dtm[0]) and (los_dtm[1] >= -los_dtm[0]):
                            # 4.2.3.2 - Vis?e principalement orient?e est
                            #         - Intersection avec le c?t? est
                            a_2 = (row_c + 1 - row_0) / los_dtm[1]
                            p_2[0] = col_0 + a_2 * los_dtm[0]

                            if (p_2[0] > col_c) and (p_2[0] < (col_c + 1)):
                                #  La vis?e sort par l'est
                                row_c += 1
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[0] < col_c:
                                # La vis?e sort par le nord
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[0] > (col_c + 1):
                                # La vis?e sort par le sud
                                col_c += 1
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[0] == col_c:
                                # La vis?e sort par le coin nord-est
                                row_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[0] == (col_c + 1):
                                # La vis?e sort par le coin sud-est
                                col_c += 1
                                row_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                        elif (los_dtm[0] > 0) and (los_dtm[0] >= los_dtm[1]) and (los_dtm[0] >= -los_dtm[1]):
                            # 4.2.3.3 - Vis?e principalement orient?e sud
                            #         - Intersection avec le c?t? sud
                            a_2 = (col_c + 1 - col_0) / los_dtm[0]
                            p_2[1] = row_0 + a_2 * los_dtm[1]

                            if (p_2[1] > row_c) and (p_2[1] < (row_c + 1)):
                                # La vis?e sort par le sud
                                col_c += 1
                                p_2[0] = col_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[1] < row_c:
                                # La vis?e sort par l'ouest
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[1] > row_c + 1:
                                # La vis?e sort par l'est
                                row_c += 1
                                a_2 = (row_c - row_0) / los_dtm[1]
                                p_2[0] = col_0 + a_2 * los_dtm[0]
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[1] == row_c:
                                # La vis?e sort par le coin sud-ouest
                                col_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[1] == row_c + 1:
                                # La vis?e sort par le coin sud-est
                                col_c += 1
                                row_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                        elif (los_dtm[1] < 0) and (los_dtm[1] <= los_dtm[0]) and (los_dtm[1] <= -los_dtm[0]):
                            #  4.2.3.4 - Vis?e principalement orient?e ouest
                            #          - Intersection avec le c?t? ouest
                            a_2 = (row_c - row_0) / los_dtm[1]
                            p_2[0] = col_0 + a_2 * los_dtm[0]

                            if (p_2[0] > col_c) and (p_2[0] < col_c + 1):
                                # La vis?e sort par l'ouest
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                            elif p_2[0] < col_c:
                                # La vis?e sort par le nord
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1

                            elif p_2[0] > (col_c + 1):
                                # La vis?e sort par le sud
                                col_c += 1
                                p_2[0] = col_c
                                a_2 = (col_c - col_0) / los_dtm[0]
                                p_2[1] = row_0 + a_2 * los_dtm[1]
                                p_2[2] = z_0 + a_2 * los_dtm[2]

                            elif p_2[0] == col_c:
                                # La vis?e sort par le coin nord-ouest
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                col_c -= 1
                                row_c -= 1

                            elif p_2[0] == (col_c + 1):
                                # La vis?e sort par le coin sud-ouest
                                col_c += 1
                                p_2[0] = col_c
                                p_2[1] = row_c
                                p_2[2] = z_0 + a_2 * los_dtm[2]
                                row_c -= 1

                    # LDD - Verification des bornes min et max de la "couche"
                    b_intersect = False

                    if p_2[2] > z_0:
                        # On est remonte trop haut, et ca c'est pas bon !!!
                        b_intersect = not (((p_1[2] > h_s) and (z_0 > h_s)) or ((p_1[2] < h_i) and (z_0 < h_i)))

                    elif p_2[2] < z_1:
                        # On est descendu trop bas, et ca c'est pas bon non plus !!! (m?me si c'est d?j? plus logique)
                        b_intersect = not (((p_1[2] > h_s) and (z_1 > h_s)) or ((p_1[2] < h_i) and (z_1 < h_i)))

                    else:
                        b_intersect = not (((p_1[2] > h_s) and (p_2[2] > h_s)) or ((p_1[2] < h_i) and (p_2[2] < h_i)))

                    # 5. Test d'intersection de la vis?e avec le cube
                    if b_intersect:
                        # Il y a intersection entre la vis?e et le cube
                        # 5.1 - Altitudes du DTM
                        alti_1 = self.interpolate(p_1[0], p_1[1])
                        h_2 = self.interpolate(p_2[0], p_2[1])

                        # 5.2 - Diff?rences d'altitude avec le DTM
                        d_alti_1 = p_1[2] - alti_1
                        d_2 = p_2[2] - h_2

                        # 5.3 - Test d'intersection avec le DTM
                        if d_alti_1 * d_2 <= 0:
                            # Il y a intersection entre la vis?e et le DTM
                            # 5.3.1 - Calcul de la solution approch?e
                            d_2 = 2 * self.tol_z  # Init de d_2 > TOL_Z
                            col_a = p_2[0]
                            row_a = p_2[1]
                            z_a = h_2

                            while abs(d_2) > self.tol_z:
                                # 5.3.1.1 - Coefficient d'interpolation lin?aire de h
                                c_h = (p_1[2] - alti_1) / ((h_2 - alti_1) - (p_2[2] - p_1[2]))

                                # 5.3.1.2 - Position du point interpole
                                col_a = p_1[0] + c_h * (p_2[0] - p_1[0])
                                row_a = p_1[1] + c_h * (p_2[1] - p_1[1])
                                z_a = p_1[2] + c_h * (p_2[2] - p_1[2])

                                # 5.3.1.3 - Altitude du point interpole
                                z_v = self.interpolate(col_a, row_a)

                                # 5.3.1.4 - Ecart d'altitude au point interpole
                                d_2 = z_v - z_a

                                # 5.3.1.5 - Mise a jour
                                if d_2 < 0:
                                    # Mise a jour du point haut
                                    p_1[0] = col_a
                                    p_1[1] = row_a
                                    p_1[2] = z_a
                                    alti_1 = z_v

                                else:
                                    # Mise ? jour du point bas
                                    p_2[0] = col_a
                                    p_2[1] = row_a
                                    p_2[2] = z_a
                                    h_2 = z_v

                            # // Fin, retour
                            p_1[0] = col_a
                            p_1[1] = row_a
                            p_1[2] = z_a

                            b_trouve = True
                            point_r = self.index_to_ter(p_1)
                            return True, b_trouve, point_r

                # Fin boucle sur les mailles

                # Test si on est toujours dans le cube DTM
                if a_2 >= 1:
                    # Changement de plan
                    i_0 += 1

                    # Chargement dans p_2 du nouveau sommet
                    p_2 = los_index[i_0].copy()
                    h_intersect_p2 = alti[i_0].copy()

                else:

                    # LDD - On a boucl? sur les mailles, on n'a rien trouv? et on n'a pas atteint le plan suivant
                    # Ca veut dire qu'on sort de l'emprise, pas la peine de continuer
                    b_trouve = False
                    return True, b_trouve, point_r
            # Fin cas general (visee non verticale)

        # Fin boucle sur les sommets
        # Fin, retour
        b_trouve = False
        return True, point_r
