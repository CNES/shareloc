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
import os

#------------------------------------------------------------------------------
class DTM:
    def __init__(self, dtm_filename, dtm_format ='bsq'):
        self.dtm_file    = dtm_filename
        self.format     = dtm_format
        self.Z          = None
        self.Zmin       = None
        self.Zmax       = None
        self.x0         = None
        self.y0         = None
        self.px         = None
        self.py         = None
        self.nc         = None
        self.nl         = None
        self.a          = None
        self.b          = None
        self.c          = None
        self.d          = None
        self.Zmin_cell  = None
        self.Zmax_cell  = None
        self.TOL_Z      = 0.0001

        #lecture mnt
        self.charge()
        self.InitMinMax()
        self.Zmax    = self.Z.max()
        self.Zmin    = self.Z.min()

        self.a = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.b = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        self.c = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        self.d = np.array([0.0, self.nl-1.0 , 0.0, self.nc-1.0, self.Zmin, self.Zmax])

        self.plans = np.array([[1.0,0.0,0.0,0.0],\
                            [1.0,0.0,0.0,self.nl-1.0],
                            [0.0,1.0,0.0,0.0,],
                            [0.0,1.0,0.0,self.nc-1.0],
                            [0.0,0.0,1.0,self.Zmin],
                            [0.0,0.0,1.0,self.Zmax]])

    def charge(self):
        if self.format=='bsq':
            hd_babel = read_hdbabel_header(self.dtm_file)
            for key in hd_babel:
                setattr(self,key,hd_babel[key])
            self.Z = read_bsq_grid(self.dtm_file, self.nl, self.nc, self.data_type)
        else:
            print("format de mnt non reconnu")


    def eq_plan(self,i,P):
        return self.a[i]*P[0] + self.b[i]*P[1] + self.c[i]*P[2] - self.d[i]

    def TerToDTM(self,vectTer):
        #VectTer : (x, y, alt dans rep sol)
        #VectDTM de taille (n,3) ou (3)
        vectDTM = vectTer.copy()
        vectDTM[0] = (vectTer[1] - self.y0) / self.py
        vectDTM[1] = (vectTer[0] - self.x0) / self.px
        return vectDTM

    def TersToDTMs(self,vectTers):
        #VectTer : (x, y, alt dans rep sol)
        #VectDTM de taille (n,3) ou (3)
        vectDTMs = vectTers.copy()
        for i,vectTer in enumerate(vectTers):
            vectDTMs[i,:] = self.TerToDTM(vectTer)
        return vectDTMs

    def DTMToTer(self,vectDTM):
        #VectDTM : (lig, col, alt dans grille MNT)
        #VectDTM de taille (n,3) ou (3)
        vectTer = vectDTM.copy()
        vectTer[0] = self.x0 + self.px * vectDTM[1]
        vectTer[1] = self.y0 + self.py * vectDTM[0]
        return vectTer

    def MakeAlti(self,dX,dY):
        """
        Indices des noeuds de la maille incluant le point
        Remarque : Si le point est en fin de MNT, on prend l'avant dernier point du MNT comme reference
                   pour pouvoir calculer son altitude avec des coefficients d'interpolation (1,1)
                   plutot que de segmentation faulter avec des coefficients d'interpolation (0,0)
        On va meme blinder plus que ca d'une maniere generale : si les coordonnees dX, dY debordent du MNT
        On "clipe" l'indice i1 pour qu'il rentre dans l'intervalle [0, _nl -2]
        On "clipe" l'indice j1 pour qu'il rentre dans l'intervalle [0, _nc -2]

        Attention on est en convention "cellule" : [0,0] est le centre de la première cellules (pas [0.5,0.5])
        """

        if (dX < 0):
            i1 = 0
        elif (dX >= self.nl-1):
            i1 = self.nl - 2
        else:
            i1 = int(np.floor(dX))

        i2 = i1+1

        if (dY < 0):
            j1 = 0
        elif (dY >= self.nc-1):
            j1 = self.nc - 2
        else:
            j1 = int(np.floor(dY))

        j2 = j1+1
        #Coefficients d'interpolation bilineaire
        u = dY - j1
        v = dX - i1
        #Altitude
        alt = (1-u)*(1-v)*self.Z[i1,j1] + u*(1-v)*self.Z[i1,j2] +\
              (1-u)*v*self.Z[i2,j1]     + u*v*self.Z[i2,j2]

        return alt

    def InitMinMax(self):
        #Generation des altitudes inf et sup de chaque maille du MNT
        NC = self.nc-1
        NL = self.nl-1

        self.Zmax_cell = np.zeros((NL,NC))
        self.Zmin_cell = np.zeros((NL,NC))

        # On calcule les altitudes min et max du MNT
        nMin =  32000
        nMax = -32000

        for i in range(NL):
            k = 0
            for j in range(NC):
                k +=1
                dAltMin = nMin
                dAltMax = nMax

                dZ1 = self.Z[i,j]
                if (dZ1 < dAltMin): dAltMin = dZ1
                if (dZ1 > dAltMax): dAltMax = dZ1

                dZ2 =  self.Z[i,j+1]
                if (dZ2 < dAltMin): dAltMin = dZ2
                if (dZ2 > dAltMax): dAltMax = dZ2

                dZ3 = self.Z[i+1,j]
                if (dZ3 < dAltMin): dAltMin = dZ3
                if (dZ3 > dAltMax): dAltMax = dZ3

                dZ4 = self.Z[i+1,j+1]
                if (dZ4 < dAltMin): dAltMin = dZ4
                if (dZ4 > dAltMax): dAltMax = dZ4

                # GDN Correction BUG Interesctor
                # Il ne faut surtout pas prendre l'arrondi pour plusieurs raisons
                # 1. l'algo ulterieur ne resiste pas touours bien lorsque la maille est plate,
                #    il est donc deconseille de fournir en sortie i_altmin = i_altmax
                #    a moins que ce soi vraiment le cas en valeurs reelles
                # 2. il ne faut pas initialiser les mailles consecutives de la meme maniere par arrondi
                #    car si l'altitude min de l'une correspond a l'altitude max de l'autre il faut les distinguer
                #    par un ceil et un floor pour que les cubes se chevauchent legerement en altitude et non pas jointifs strictement
                i_altmin = np.floor(dAltMin)
                i_altmax = np.ceil(dAltMax)
                self.Zmin_cell[i,j] = i_altmin
                self.Zmax_cell[i,j] = i_altmax
        return True