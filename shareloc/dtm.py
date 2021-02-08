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
from shareloc.readwrite import read_hdbabel_header, read_bsq_grid


class DTM:
    """ DTM class dedicated
     works only with BSQ file format
     we work in cell convention [0,0] is the first cell center (not [0.5,0.5])
    """
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

        #lecture mnt
        self.load()
        self.init_min_max()
        self.alt_max = self.alt_data.max()
        self.alt_min = self.alt_data.min()

        self.plane_coef_a = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.plane_coef_b = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        self.plane_coef_c = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        self.plane_coef_d = np.array([0.0, self.row_nb-1.0, 0.0, self.column_nb-1.0, self.alt_min, self.alt_max])

        self.plans = np.array([[1.0, 0.0, 0.0, 0.0],
                               [1.0, 0.0, 0.0, self.row_nb-1.0],
                               [0.0, 1.0, 0.0, 0.0],
                               [0.0, 1.0, 0.0, self.column_nb-1.0],
                               [0.0, 0.0, 1.0, self.alt_min],
                               [0.0, 0.0, 1.0, self.alt_max]])

    def load(self):
        """
        load DTM infos
        """
        if self.format == "bsq":
            hd_babel = read_hdbabel_header(self.dtm_file)
            for key in hd_babel:
                setattr(self, key, hd_babel[key])
            self.alt_data = read_bsq_grid(self.dtm_file, self.row_nb, self.column_nb, self.data_type)
        else:
            print("dtm format not handled")

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
        return self.plane_coef_a[i]*position[0] + self.plane_coef_b[i]*position[1] + \
               self.plane_coef_c[i]*position[2] - self.plane_coef_d[i]

    def ter_to_index(self, vect_ter):
        """
        terrain to index conversion
        :param vect_ter: terrain coordinate
        :type vect_ter: array (1x2 or 1x3) if dimension is 3 , last coordinate is unchanged (alt)
        :return index coordinates
        :rtype array (1x2 or 1x3)
        """
        vect_dtm = vect_ter.copy()
        vect_dtm[0] = (vect_ter[1] - self.origin_y) / self.pixel_size_y
        vect_dtm[1] = (vect_ter[0] - self.origin_x) / self.pixel_size_x
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
        vect_ter[0] = self.origin_x + self.pixel_size_x * vect_dtm[1]
        vect_ter[1] = self.origin_y + self.pixel_size_y * vect_dtm[0]
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

        #index clipping in [0, _row_nb -2]
        if pos_row < 0:
            row_inf = 0
        elif pos_row >= self.row_nb-1:
            row_inf = self.row_nb - 2
        else:
            row_inf = int(np.floor(pos_row))

        row_sup = row_inf + 1

        #index clipping in [0, _nc -2]
        if pos_col < 0:
            col_inf = 0
        elif pos_col >= self.column_nb-1:
            col_inf = self.column_nb - 2
        else:
            col_inf = int(np.floor(pos_col))

        col_sup = col_inf + 1
        #Coefficients d'interpolation bilineaire
        delta_u = pos_col - col_inf
        delta_v = pos_row - row_inf
        #Altitude
        alt = (1-delta_u)*(1-delta_v)*self.alt_data[row_inf, col_inf] + \
              delta_u*(1-delta_v)*self.alt_data[row_inf, col_sup] +\
              (1-delta_u)*delta_v*self.alt_data[row_sup, col_inf] + delta_u*delta_v*self.alt_data[row_sup, col_sup]
        return alt

    def init_min_max(self):
        """
        initialize min/max at each dtm cell
        """
        column_nb_minus = self.column_nb-1
        row_nb_minus = self.row_nb-1

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

                d_z2 = self.alt_data[i, j+1]
                if d_z2 < d_alt_min:
                    d_alt_min = d_z2
                if d_z2 > d_alt_max:
                    d_alt_max = d_z2

                d_z3 = self.alt_data[i+1, j]
                if d_z3 < d_alt_min:
                    d_alt_min = d_z3
                if d_z3 > d_alt_max:
                    d_alt_max = d_z3

                d_z4 = self.alt_data[i+1, j+1]
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

    def checkCubeDTM(self,LOS):
        """
        DTM cube intersection
        :param LOS :  line of sight
        :type LOS : numpy.array
        :return intersection information (True,an intersection has been found ?, (lon,lat) of dtm position, altitude)
        :rtype tuple (bool, bool, numpy.array, float)
        """
        #LOS: (n,3):
        PointB = None
        dH3D = None
        (ui, vi, zi, hi) = ([], [], [], [])
        nbalt = LOS.shape[0]
        LOSDTM = self.ters_to_indexs(LOS)
        # -----------------------------------------------------------------------
        # Nombre d'intersections valides trouvees
        nbi = 0
        # -----------------------------------------------------------------------
        # On boucle sur les plans du cube DTM
        for f in range(6):
            # -----------------------------------------------------------------------
            # Initialisation du sommet de la visee
            sB = LOSDTM[0,:]
            # -----------------------------------------------------------------------
            # Initialisation de la position par / au plan
            #print self.plans[f,:-1]
            #posB = (self.plans[f,:-1]*sB).sum() - self.d[f]
            #print posB
            posB = self.eq_plan(f,sB)
            # -----------------------------------------------------------------------
            # On boucle sur les segments de la visee et on controle
            # si on traverse ou pas la face courante f du cube DTM
            for p in range(nbalt):
                # -----------------------------------------------------------------------
                # Transfert du point B dans le point A
                posA = posB
                sA = sB.copy()
                # -----------------------------------------------------------------------
                # Reinit du point B
                sB = LOSDTM[p,:] #dg on itere sur les differents points de la visee
                #print sB
                # -----------------------------------------------------------------------
                # Initialisation de la position par / au plan
                posB = self.eq_plan(f,sB)
                #print 'posAposB',posA,posB
                # -----------------------------------------------------------------------
                # Test d'intersection : posA et posB de signes opposes
                if posA * posB <= 0:
                    # -----------------------------------------------------------------------
                    if not posA and not posB:
                        # Trop de solutions !! (A et B sont sur le plan) #comment on le gere ??
                        continue
                    # -----------------------------------------------------------------------
                    elif not posA :
                        # A est solution (il est sur le plan)
                        ui.append(sA[0])
                        vi.append(sA[1])
                        zi.append(sA[2])
                        hi.append(p - 1)
                    # -----------------------------------------------------------------------
                    elif not posB:
                        # B est solution (il est sur le plan)
                        ui.append(sB[0])
                        vi.append(sB[1])
                        zi.append(sB[2])
                        hi.append(p)
                    # -----------------------------------------------------------------------
                    else:
                        # -----------------------------------------------------------------------
                        # A et B sont de part et d'autre du plan
                        # Coefficients d'interpolation de l'intersection
                        # entre A et B
                        cA = posB / (posB - posA)
                        cB = -posA / (posB - posA)
                        # Affectation ou interpolation
                        # NB : pour eviter les pb lors du test
                        #      <estSurCube> (voir + loin)
                        # . coordonn?e <u> (ligne)
                        # -----------------------------------------------------------------------
                        if f < 2:
                            ui.append(self.plane_coef_d[f])
                        # -----------------------------------------------------------------------
                        else:
                            ui.append(cA * sA[0] + cB * sB[0])
                        # -----------------------------------------------------------------------
                        # . coordonnee <v> (colonne)
                        if f > 1 and f < 4:
                            vi.append(self.plane_coef_d[f])
                        else:
                            vi.append(cA * sA[1] + cB * sB[1])
                        # -----------------------------------------------------------------------
                        # . coordonn?e <z> (altitude)
                        if f > 3:
                            zi.append(self.plane_coef_d[f])
                        # -----------------------------------------------------------------------
                        else:
                            zi.append(cA * sA[2] + cB * sB[2])
                        # . coordonn?e <h> (abscisse visee)
                        hi.append(p - cA) #index non entier de l'intersection
                    # -----------------------------------------------------------------------
                    # Incrementation du nombre d'intersections trouvees
                    nbi += 1
                    # -----------------------------------------------------------------------
                    # Passage a la face du cube suivante
                    break

        # -----------------------------------------------------------------------
        # Tri des points le long de la visee (il y en a au moins deux)
        #on les range par ordre decroissant en altitude
        for p in range(nbi):
            for q in range(p+1,nbi):
                if hi[q] < hi[p]:
                    dtmp = ui[p]
                    ui[p] = ui[q]
                    ui[q] = dtmp
                    dtmp = vi[p]
                    vi[p] = vi[q]
                    vi[q] = dtmp
                    dtmp = zi[p]
                    zi[p] = zi[q]
                    zi[q] = dtmp
                    dtmp = hi[p]
                    hi[p] = hi[q]
                    hi[q] = dtmp

        # -----------------------------------------------------------------------
        # Filtrage des points non situes sur le cube
        p = 0
        while p < nbi:
            #test a l'interieur du cube
            test_on_cube = (ui[p] >= self.plane_coef_d[0]) and (ui[p] <= self.plane_coef_d[1]) and \
                         (vi[p] >= self.plane_coef_d[2]) and (vi[p] <= self.plane_coef_d[3]) and \
                         (zi[p] >= self.plane_coef_d[4]) and (zi[p] <= self.plane_coef_d[5])
            if not test_on_cube:
                # On translate tous les points suivants (on ecrase ce point non valide)
                for q in range(p+1, nbi):
                    ui[q - 1] = ui[q]
                    vi[q - 1] = vi[q]
                    zi[q - 1] = zi[q]
                    hi[q - 1] = hi[q]
                nbi -= 1
            else:
                p += 1
        # -----------------------------------------------------------------------
        # Pas de solution si 0 ou 1 seul point trouve (on a tangente le cube)
        if nbi < 2:
            bTrouve = False
            return (True, bTrouve, PointB, dH3D)
        # -----------------------------------------------------------------------
        # Il ne reste que 2 points donc on traverse le cube
        # LAIG-FA-MAJA-2168-CNES: plus de filtrage sur les point identiques. Il peut y avoir un nombre depoints > 2
        # Initialisation du point courant
        # Coordonnees DTM
        point_dtm = np.zeros(3)
        point_dtm[0] = ui[0]
        point_dtm[1] = vi[0]
        point_dtm[2] = zi[0]
        #PointDTM est la premiere intersection avec le cube (lig, col)
        # -----------------------------------------------------------------------
        # h dans gld 3D
        dH3D = hi[0]
        #dH3D correspond a l'index (non entier) d'interpolation en h
        # -----------------------------------------------------------------------
        # Coordonnees terrain
        point_b = self.index_to_ter(point_dtm)
        #PointB est le point Terrain (lon,lat)
        # -----------------------------------------------------------------------
        # Fin, retour
        bTrouve = True
        return (True, bTrouve, point_b,dH3D)

    def intersection(self, LOS, PointB, dH3D):
        """
        DTM intersection
        :param LOS :  line of sight
        :type LOS : numpy.array
        :param PointB :  position of intersection in DTM cube
        :type PointB : numpy.array
        :param dH3D :  altitude in DTM cube
        :type dH3D : float
        :return intersection information (True,an intersection has been found ?, position of intersection)
        :rtype tuple (bool, bool, numpy.array)
        """
        LOSDTM = self.ters_to_indexs(LOS)
        PointB_DTM = self.ter_to_index(PointB)
        PointR = np.zeros(3)
        (npl, _) = LOS.shape
        H3D = range(npl, -1, -1)

        p1 = PointB_DTM.copy() #[p1[0],p1[1],p1[2]]

        dH3D_p1 = dH3D

        nu = self.row_nb
        nv = self.column_nb

        #1 - Initilialisation et tests prealables
        #   1.1 - Test si le sommet est au-dessu du DTM
        #       - Calcul de l'altitude du DTM ? la position du sommet
        h1 = self.interpolate(p1[0], p1[1])
        #       - Calcul de l'ecart d'altitude au DTM
        d1 = p1[2] - h1

        #       - Test si le nouveau point haut est au dessus du DTM
        if d1 < 0:
        #       - Point situ? en dessous du DTM
        #          . ceci signifie que la vis?e rentre dans le DTM par le c?t?
        #          . donc en dessous, pas de solution
            bTrouve = False
            return (True, bTrouve, PointR)


        #   1.2 - Initialisation du rang du premier sommet de la visee
        i0 = int(np.floor(dH3D_p1))

        #   1.3 - Initialisation du point de depart (dans p2)
        p2 = PointB_DTM.copy()
        dH3D_p2 = dH3D

        #2. - Boucle sur les plans de grille
        while i0 < (LOSDTM.size - 1):
            #2.1 - Initialisation du sommet courant de la visee
            u0 = LOSDTM[i0][0]
            v0 = LOSDTM[i0][1]
            z0 = LOSDTM[i0][2]
            z1 = LOSDTM[i0 + 1][2]

            #2.2 - Initialisation de la visee DTM
            VM = LOSDTM[i0 + 1] - LOSDTM[i0]

            #2.3 - Test si visee verticale
            if (VM[0]==0 and VM[1]==0):
                #2.3.1 - La vis?e est verticale :
                #    - Calcul de l'altitude du DTM ? la position du sommet
                h1 = self.interpolate(u0, v0)

                #    Test si le plan suivant est en dessous du DTM
                if LOSDTM[i0 + 1][2] <= h1:
                    #Init point de sortie
                    p1[0] = u0
                    p1[1] = v0
                    p1[2] = h1
                    bTrouve = True
                    self.index_to_ter(p1, PointR)
                    return (True, bTrouve, PointR)
                else:
                    #Positionnement sur le sommet suivant
                    i0 += 1
            else:
                #2.3.2 - La visee n'est pas verticale :
                #         elle va donc survoler le DTM
                #         . on peut donc poursuivre
                #         . reste ? d?montrer que la vis?e se crashe sur le DTM...
                #
                # Initialisation du point de d?part
                # Initialisation de son abscisse sur la vis?e
                a2 = dH3D_p2 - i0

                # Correction d'un bug FA DG 10, a la reinitialisation a2 peut valoir 1
                # Ce qui a pour consequence de sauter au segment suivant de la visee
                # Ainsi, on peut perdre des points sur une tranche d'altitude
                # Pour etre sur de scanner le segment de visee, on reinitialise
                # le point au depart du segment, soit a2 = 0
                if a2 >= 1.:
                    a2 = 0.

                # Initialisation de la premi?re maille DTM intersect?e
                #  - Initialisation des indices de la maille
                uc = int(np.floor(p2[0]))
                vc = int(np.floor(p2[1]))

                # NB :    pr?caution avant de d?marrer :
                #        . on se met du bon cote de la maille
                #        . en principe, on ne doit pas sortir du DTM
                # On rentre par le bas, la maille DTM est la precedente
                if (p2[0] == uc) and (VM[0] < 0):
                    uc -= 1

                # On rentre par la gauche, la maille DTM est la precedente
                if (p2[1] == vc) and (VM[1] < 0):
                    vc -= 1

                # LDD - On est deja en dehors des limites, on s'arrete
                if not((a2 < 1) and (uc > -1) and (uc < (nu - 1)) and (vc > -1) and (vc < (nv - 1)) ) and (a2 < 1):
                    bTrouve = False
                    (True, bTrouve, PointR)

                # Boucle de recherche iterative de la maille intersectee
                while (a2 < 1) and (uc > -1) and (uc < (nu - 1)) and (vc > -1) and (vc < (nv - 1)):
                    # - Altitudes min et max de la maille
                    hi = self.alt_min_cell[uc, vc]
                    hs = self.alt_max_cell[uc, vc]

                    # - Transfert : le point bas devient le point haut
                    #a1 = a2;
                    # p1 devient p2
                    p1 = p2.copy()
                    dH3D_p1 = dH3D_p2.copy()

                    # 4.2 - D?termination d'un nouveau point bas
                    #      - Test d'orientation de la vis?e
                    if not VM[0]:
                        # 4.2.1 - La vis?e est orientee pile poil est-ouest
                        #   p2[0] = p1[0] ; // inutile, est deja initialise
                        #       - Test d'orientation de la visee
                        if VM[1] < 0:
                            # 4.2.1.1 - La vis?e part plein ouest
                            p2[1] = vc
                            vc -= 1
                        else:
                            # 4.2.1.2 - La vis?e part plein est
                            vc += 1
                            p2[1] = vc

                        a2 = (p2[1] - v0) / VM[1]
                        p2[2] = z0 + a2 * VM[2]

                    elif not VM[1]:
                        # 4.2.2 - La vis?e est orient?e nord-sud
                        #  p2[1] = p1[1] ;
                        #       - Test d'orientation de la visee
                        if VM[0] < 0:
                            # 4.2.2.1 - La vis?e part plein nord
                            p2[0] = uc
                            uc -= 1
                        else:
                            # 4.2.2.2 - La vis?e part plein sud
                            uc += 1
                            p2[0] = uc

                        a2 = (p2[0] - u0) / VM[0]
                        p2[2] = z0 + a2 * VM[2]
                    else:
                        # 4.2.3 - La vis?e est quelconque
                        #            - D?termination du cot? de sortie
                        if (VM[0] < 0) and (VM[0] <= VM[1]) and (VM[0] <= -VM[1]):
                            # 4.2.3.1 - Vis?e principalement orient?e nord
                            #             - Intersection avec le c?t? nord
                            a2    = (uc - u0) / VM[0]
                            p2[1] = v0 + a2 * VM[1]

                            if (p2[1] > vc) and (p2[1] < (vc + 1)):
                                # La vis?e sort par le nord
                                p2[0] = uc
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1

                            elif p2[1] < vc:
                                # La vis?e sort par l'ouest
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc -= 1

                            elif p2[1] > (vc + 1):
                                # La vis?e sort par l'est
                                vc += 1
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                            elif p2[1] == vc:
                                # La vis?e sort par le coin nord-ouest
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc -= 1
                                vc -= 1

                            elif p2[1] == (vc + 1):
                                # La vis?e sort par le coin nord-est
                                p2[0] = uc
                                vc += 1
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc -= 1

                        elif (VM[1] > 0) and (VM[1] >= VM[0]) and (VM[1] >= -VM[0]):
                            # 4.2.3.2 - Vis?e principalement orient?e est
                            #         - Intersection avec le c?t? est
                            a2 = (vc + 1 - v0) / VM[1]
                            p2[0] = u0 + a2 * VM[0]

                            if (p2[0] > uc) and (p2[0] < (uc + 1)):
                                #  La vis?e sort par l'est
                                vc+=1
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                            elif p2[0] < uc:
                                # La vis?e sort par le nord
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]
                                uc -= 1

                            elif p2[0] > (uc + 1):

                                # La vis?e sort par le sud
                                uc += 1
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]

                            elif p2[0] == uc:
                                # La vis?e sort par le coin nord-est
                                vc += 1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc -= 1

                            elif p2[0] == (uc + 1):
                                # La vis?e sort par le coin sud-est
                                uc += 1
                                vc += 1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                        elif (VM[0] > 0) and (VM[0] >= VM[1]) and (VM[0] >= -VM[1]):
                            # 4.2.3.3 - Vis?e principalement orient?e sud
                            #         - Intersection avec le c?t? sud
                            a2 = (uc + 1 - u0) / VM[0]
                            p2[1] = v0 + a2 * VM[1]

                            if (p2[1] > vc) and (p2[1] < (vc + 1)):
                                # La vis?e sort par le sud
                                uc += 1
                                p2[0] = uc
                                p2[2] = z0 + a2 * VM[2]

                            elif p2[1] < vc:
                                # La vis?e sort par l'ouest
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc -= 1

                            elif p2[1] > vc + 1:
                                # La vis?e sort par l'est
                                vc+=1
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                            elif p2[1] == vc:
                                # La vis?e sort par le coin sud-ouest
                                uc += 1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc -= 1

                            elif p2[1] == vc + 1:
                                # La vis?e sort par le coin sud-est
                                uc += 1
                                vc += 1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                        elif (VM[1] < 0) and (VM[1] <= VM[0]) and (VM[1] <= -VM[0]):
                            #  4.2.3.4 - Vis?e principalement orient?e ouest
                            #          - Intersection avec le c?t? ouest
                            a2 = (vc - v0) / VM[1]
                            p2[0] = u0 + a2 * VM[0]

                            if (p2[0] > uc) and (p2[0] < uc + 1):
                                # La vis?e sort par l'ouest
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc -= 1

                            elif p2[0] < uc:
                                # La vis?e sort par le nord
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]
                                uc -= 1

                            elif p2[0] > (uc + 1):
                                # La vis?e sort par le sud
                                uc += 1
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]

                            elif p2[0] == uc:
                                # La vis?e sort par le coin nord-ouest
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc -= 1
                                vc -= 1

                            elif p2[0] == (uc + 1):
                                # La vis?e sort par le coin sud-ouest
                                uc += 1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc -= 1

                    # LDD - Verification des bornes min et max de la "couche"
                    bIntersect = False

                    if p2[2] > z0:
                        # On est remonte trop haut, et ca c'est pas bon !!!
                        bIntersect = not(((p1[2] > hs) and (z0 > hs)) or ((p1[2] < hi) and (z0 < hi)))

                    elif p2[2] < z1:
                        # On est descendu trop bas, et ca c'est pas bon non plus !!! (m?me si c'est d?j? plus logique)
                        bIntersect = not(((p1[2] > hs) and (z1 > hs)) or ((p1[2] < hi) and (z1 < hi)))

                    else:
                        bIntersect = not(((p1[2] > hs) and (p2[2] > hs)) or ((p1[2] < hi) and (p2[2] < hi)))

                    # 5. Test d'intersection de la vis?e avec le cube
                    if bIntersect:
                        # Il y a intersection entre la vis?e et le cube
                        # 5.1 - Altitudes du DTM
                        h1 = self.interpolate(p1[0], p1[1])
                        h2 = self.interpolate(p2[0], p2[1])

                        # 5.2 - Diff?rences d'altitude avec le DTM
                        d1 = p1[2] - h1
                        d2 = p2[2] - h2

                        # 5.3 - Test d'intersection avec le DTM
                        if d1 * d2 <= 0:
                            # Il y a intersection entre la vis?e et le DTM
                            # 5.3.1 - Calcul de la solution approch?e
                            d2 = 2 * self.tol_z # Init de d2 > TOL_Z
                            ua = p2[0]
                            va = p2[1]
                            za = h2

                            while abs(d2) > self.tol_z:
                                # 5.3.1.1 - Coefficient d'interpolation lin?aire de h
                                ch = (p1[2] - h1) / ((h2 - h1) - (p2[2] - p1[2]))

                                # 5.3.1.2 - Position du point interpole
                                ua = p1[0] + ch * (p2[0] - p1[0])
                                va = p1[1] + ch * (p2[1] - p1[1])
                                za = p1[2] + ch * (p2[2] - p1[2])

                                # 5.3.1.3 - Altitude du point interpole
                                zv = self.interpolate(ua, va)

                                # 5.3.1.4 - Ecart d'altitude au point interpole
                                d2 = zv - za

                                # 5.3.1.5 - Mise a jour
                                if d2 < 0:
                                    # Mise a jour du point haut
                                    p1[0] = ua
                                    p1[1] = va
                                    p1[2] = za
                                    h1 = zv

                                else:
                                    # Mise ? jour du point bas
                                    p2[0] = ua
                                    p2[1] = va
                                    p2[2] = za
                                    h2 = zv

                            #// Fin, retour
                            p1[0] = ua
                            p1[1] = va
                            p1[2] = za

                            bTrouve = True
                            PointR = self.index_to_ter(p1)
                            return (True, bTrouve, PointR)

                # Fin boucle sur les mailles

                # Test si on est toujours dans le cube DTM
                if a2 >= 1:
                    # Changement de plan
                    i0 += 1

                    # Chargement dans p2 du nouveau sommet
                    p2 = LOSDTM[i0].copy()
                    dH3D_p2 = H3D[i0].copy()

                else:

                    # LDD - On a boucl? sur les mailles, on n'a rien trouv? et on n'a pas atteint le plan suivant
                    # Ca veut dire qu'on sort de l'emprise, pas la peine de continuer
                    bTrouve = False
                    return (True, bTrouve, PointR)
            # Fin cas general (visee non verticale)

        # Fin boucle sur les sommets

        # Fin, retour
        bTrouve = False
        return (True, PointR)
