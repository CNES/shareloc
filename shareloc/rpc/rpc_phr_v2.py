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

from xml.dom import minidom
from numpy import array, dot, zeros, sqrt
from os.path import basename

def renvoie_linesep(txt_liste_lines):
	"""Renvoie le separateur de ligne d'un texte sous forme de liste de lignes
	Obtenu par readlines
	"""
	if txt_liste_lines[0].endswith('\r\n'):
		line_sep = '\r\n'
	elif txt_liste_lines[0].endswith('\n'):
		line_sep = '\n'
	return line_sep


def read_eucl_file(eucl_file):
    """
    read euclidium file and parse it
    :param eucl_file : euclidium file
    :type eucl_file : str
    :return parsed file
    :rtype dict
    """
    parsed_file = dict()
    with open(eucl_file, 'r') as fid:
        txt = fid.readlines()

    lsep = renvoie_linesep(txt)

    ind_debut_PX = txt.index('>>\tCOEFF POLYNOME PXOUT' + lsep)
    ind_debut_QX = txt.index('>>\tCOEFF POLYNOME QXOUT' + lsep)
    ind_debut_PY = txt.index('>>\tCOEFF POLYNOME PYOUT' + lsep)
    ind_debut_QY = txt.index('>>\tCOEFF POLYNOME QYOUT' + lsep)

    coeff_PX_str = txt[ind_debut_PX + 1:ind_debut_PX + 21]
    coeff_QX_str = txt[ind_debut_QX + 1:ind_debut_QX + 21]
    coeff_PY_str = txt[ind_debut_PY + 1:ind_debut_PY + 21]
    coeff_QY_str = txt[ind_debut_QY + 1:ind_debut_QY + 21]

    parsed_file['coeff_PX'] = [float(coeff.split()[1]) for coeff in coeff_PX_str]
    parsed_file['coeff_QX'] = [float(coeff.split()[1]) for coeff in coeff_QX_str]
    parsed_file['coeff_PY'] = [float(coeff.split()[1]) for coeff in coeff_PY_str]
    parsed_file['coeff_QY'] = [float(coeff.split()[1]) for coeff in coeff_QY_str]

    # list [offset , scale]
    normalisation_coeff = dict()
    for l in txt:
        if l.startswith('>>\tTYPE_OBJET'):
            if l.split()[-1].endswith('Inverse'):
                parsed_file['type_fic'] = 'I'
            if l.split()[-1].endswith('Directe'):
                parsed_file['type_fic'] = 'D'
        if l.startswith('>>\tXIN_OFFSET'):
            lsplit = l.split()
            normalisation_coeff['X'] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tYIN_OFFSET'):
            lsplit = l.split()
            normalisation_coeff['Y'] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tZIN_OFFSET'):
            lsplit = l.split()
            normalisation_coeff['ALT'] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tXOUT_OFFSET'):
            lsplit = l.split()
            normalisation_coeff['COL'] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tYOUT_OFFSET'):
            lsplit = l.split()
            normalisation_coeff['LIG'] = [float(lsplit[4]), float(lsplit[5])]
    parsed_file['normalisation_coeffs'] = normalisation_coeff
    return parsed_file


def check_coeff_consistency(dict1, dict2):
    for key, value in dict1.items() :
        if dict2[key] != value:
            print("normalisation coeffs are different between"
                  " direct en inverse one : {} : {} {}".format(key,value,dict2[key]))


class FonctRatD:
    def __init__(self,rpc_params):
        for a, b in rpc_params.items():
            setattr(self, a, b)

        self.type = 'rpc'
        self.lim_extrapol = 1.0001
        #chaque mononome: c[0]*X**c[1]*Y**c[2]*Z**c[3]
        ordre_monomes_LAI = \
                [[1,0,0,0],[1,1,0,0],[1,0,1,0],\
                 [1,0,0,1],[1,1,1,0],[1,1,0,1],\
                 [1,0,1,1],[1,2,0,0],[1,0,2,0],\
                 [1,0,0,2],[1,1,1,1],[1,3,0,0],\
                 [1,1,2,0],[1,1,0,2],[1,2,1,0],\
                 [1,0,3,0],[1,0,1,2],[1,2,0,1],\
                 [1,0,2,1],[1,0,0,3]]

        self.Monomes    = ordre_monomes_LAI

        #coefficient des degres monomes avec derivation 1ere variable
        self.monomes_deriv_1 = \
                [[0,0,0,0],[1,0,0,0],[0,0,1,0],\
                 [0,0,0,1],[1,0,1,0],[1,0,0,1],\
                 [0,0,1,1],[2,1,0,0],[0,0,2,0],\
                 [0,0,0,2],[1,0,1,1],[3,2,0,0],\
                 [1,0,2,0],[1,0,0,2],[2,1,1,0],\
                 [0,0,3,0],[0,0,1,2],[2,1,0,1],\
                 [0,0,2,1],[0,0,0,3]]

        #coefficient des degres monomes avec derivation 1ere variable
        self.monomes_deriv_2 = \
                [[0,0,0,0],[0,1,0,0],[1,0,0,0],\
                 [0,0,0,1],[1,1,0,0],[0,1,0,1],\
                 [1,0,0,1],[0,2,0,0],[2,0,1,0],\
                 [0,0,0,2],[1,1,0,1],[0,3,0,0],\
                 [2,1,1,0],[0,1,0,2],[1,2,0,0],\
                 [3,0,2,0],[1,0,0,2],[0,2,0,1],\
                 [2,0,1,1],[0,0,0,3]]

    @classmethod
    def from_dimap(cls, dimap_filepath):
        """ load from dimap """

        rpc_params = dict()
        rpc_params['driver_type'] = 'dimapV1'
        if not basename(dimap_filepath).endswith('XML'.upper()):
            raise ValueError("dimap must ends with .xml")

        xmldoc= minidom.parse(dimap_filepath)
        #TODO check dimap version

        GLOBAL_RFM    = xmldoc.getElementsByTagName('Global_RFM')
        RFM_Validity     = xmldoc.getElementsByTagName('RFM_Validity')
        coeff_LON = [float(el) for el in GLOBAL_RFM[0].getElementsByTagName('F_LON')[0].firstChild.data.split()]
        coeff_LAT = [float(el) for el in GLOBAL_RFM[0].getElementsByTagName('F_LAT')[0].firstChild.data.split()]
        coeff_COL = [float(el) for el in GLOBAL_RFM[0].getElementsByTagName('F_COL')[0].firstChild.data.split()]
        coeff_LIG = [float(el) for el in GLOBAL_RFM[0].getElementsByTagName('F_ROW')[0].firstChild.data.split()]

        A_lon = float(RFM_Validity[0].getElementsByTagName('Lon')[0].getElementsByTagName('A')[0].firstChild.data)
        B_lon = float(RFM_Validity[0].getElementsByTagName('Lon')[0].getElementsByTagName('B')[0].firstChild.data)
        A_lat = float(RFM_Validity[0].getElementsByTagName('Lat')[0].getElementsByTagName('A')[0].firstChild.data)
        B_lat = float(RFM_Validity[0].getElementsByTagName('Lat')[0].getElementsByTagName('B')[0].firstChild.data)
        A_alt = float(RFM_Validity[0].getElementsByTagName('Alt')[0].getElementsByTagName('A')[0].firstChild.data)
        B_alt = float(RFM_Validity[0].getElementsByTagName('Alt')[0].getElementsByTagName('B')[0].firstChild.data)
        A_col = float(RFM_Validity[0].getElementsByTagName('Col')[0].getElementsByTagName('A')[0].firstChild.data)
        B_col = float(RFM_Validity[0].getElementsByTagName('Col')[0].getElementsByTagName('B')[0].firstChild.data)
        A_row = float(RFM_Validity[0].getElementsByTagName('Row')[0].getElementsByTagName('A')[0].firstChild.data)
        B_row = float(RFM_Validity[0].getElementsByTagName('Row')[0].getElementsByTagName('B')[0].firstChild.data)


        rpc_params['offset_COL']    = B_col
        rpc_params['scale_COL']    = A_col
        rpc_params['offset_LIG']    = B_row
        rpc_params['scale_LIG']    = A_row
        rpc_params['offset_ALT']    = B_alt
        rpc_params['scale_ALT']    = A_alt
        rpc_params['offset_X']    = B_lon
        rpc_params['scale_X']    = A_lon
        rpc_params['offset_Y']    = B_lat
        rpc_params['scale_Y']    = A_lat
        rpc_params['Num_X']    = coeff_LON[0:20]
        rpc_params['Den_X']    = coeff_LON[20::]
        rpc_params['Num_Y']    = coeff_LAT[0:20]
        rpc_params['Den_Y']    = coeff_LAT[20::]
        rpc_params['Num_COL']    = coeff_COL[0:20]
        rpc_params['Den_COL']    = coeff_COL[20::]
        rpc_params['Num_LIG']    = coeff_LIG[0:20]
        rpc_params['Den_LIG']    = coeff_LIG[20::]
        return cls(rpc_params)

    @classmethod
    def from_euclidium(cls, inverse_euclidium_coeff, direct_euclidium_coeff=None):
        """ load from euclidium """

        rpc_params = dict()
        rpc_params['driver_type'] = 'euclidium'

        #lecture fichier euclide
        inverse_coeffs = read_eucl_file(inverse_euclidium_coeff)

        if inverse_coeffs['type_fic'] != 'I':
            print("inverse euclidium file is of {} type".format(inverse_coeffs['type_fic']))

        rpc_params['Num_COL'] = inverse_coeffs['coeff_PX']
        rpc_params['Den_COL'] = inverse_coeffs['coeff_QX']
        rpc_params['Num_LIG'] = inverse_coeffs['coeff_PY']
        rpc_params['Den_LIG'] = inverse_coeffs['coeff_QY']

        rpc_params['normalisation_coeffs'] = inverse_coeffs['normalisation_coeffs']
        for key, value in inverse_coeffs['normalisation_coeffs'].items():
            rpc_params['offset_' + key] = value[0]
            rpc_params['scale_' + key] = value[1]

        if direct_euclidium_coeff is not None :
            direct_coeffs = read_eucl_file(direct_euclidium_coeff)
            if direct_coeffs['type_fic'] != 'D':
                print("direct euclidium file is of {} type".format(direct_coeffs['type_fic']))

            check_coeff_consistency(inverse_coeffs['normalisation_coeffs'], direct_coeffs['normalisation_coeffs'])
            rpc_params['Num_X'] = direct_coeffs['coeff_PX']
            rpc_params['Den_X'] = direct_coeffs['coeff_QX']
            rpc_params['Num_Y'] = direct_coeffs['coeff_PY']
            rpc_params['Den_Y'] = direct_coeffs['coeff_QY']
        else:
            rpc_params['Num_X'] = None
            rpc_params['Den_X'] = None
            rpc_params['Num_Y'] = None
            rpc_params['Den_Y'] = None
        print(rpc_params)
        return cls(rpc_params)





    def calcule_derivees_inv(self,lon,lat,alt):
        """ calcul analytiques des derivees partielles de la loc inverse
            DCdx: derivee de loc_inv_C p/r a X
            DLdy: derivee de loc_inv_L p/r a Y
        """

        if self.Num_COL:
            Xnorm = (lon - self.offset_X)/self.scale_X
            Ynorm = (lat - self.offset_Y)/self.scale_Y
            Znorm = (alt - self.offset_ALT)/self.scale_ALT
            monomes = array([self.Monomes[i][0]*\
                 Xnorm**int(self.Monomes[i][1])*\
                 Ynorm**int(self.Monomes[i][2])*\
                 Znorm**int(self.Monomes[i][3]) for i in range(self.Monomes.__len__())])
            NumDC = dot(array(self.Num_COL),monomes)
            DenDC = dot(array(self.Den_COL),monomes)
            NumDL = dot(array(self.Num_LIG),monomes)
            DenDL = dot(array(self.Den_LIG),monomes)

            monomes_deriv_x = array([self.monomes_deriv_1[i][0]*\
                Xnorm**int(self.monomes_deriv_1[i][1])*\
                Ynorm**int(self.monomes_deriv_1[i][2])*\
                Znorm**int(self.monomes_deriv_1[i][3]) for i in range(self.monomes_deriv_1.__len__())])

            monomes_deriv_y = array([self.monomes_deriv_2[i][0]*\
                Xnorm**int(self.monomes_deriv_2[i][1])*\
                Ynorm**int(self.monomes_deriv_2[i][2])*\
                Znorm**int(self.monomes_deriv_2[i][3]) for i in range(self.monomes_deriv_2.__len__())])

            NumDCdx = dot(array(self.Num_COL),monomes_deriv_x)
            DenDCdx = dot(array(self.Den_COL),monomes_deriv_x)
            NumDLdx = dot(array(self.Num_LIG),monomes_deriv_x)
            DenDLdx = dot(array(self.Den_LIG),monomes_deriv_x)

            NumDCdy = dot(array(self.Num_COL),monomes_deriv_y)
            DenDCdy = dot(array(self.Den_COL),monomes_deriv_y)
            NumDLdy = dot(array(self.Num_LIG),monomes_deriv_y)
            DenDLdy = dot(array(self.Den_LIG),monomes_deriv_y)

            #derive (u/v)' = (u'v - v'u)/(v*v)
            DCdx = self.scale_COL/self.scale_X*(NumDCdx*DenDC - DenDCdx*NumDC)/DenDC**2
            DCdy = self.scale_COL/self.scale_Y*(NumDCdy*DenDC - DenDCdy*NumDC)/DenDC**2
            DLdx = self.scale_LIG/self.scale_X*(NumDLdx*DenDL - DenDLdx*NumDL)/DenDL**2
            DLdy = self.scale_LIG/self.scale_Y*(NumDLdy*DenDL - DenDLdy*NumDL)/DenDL**2

        return (DCdx,DCdy,DLdx,DLdy)


    def direct_loc_dtm(self, row, col, dtm):
        """
        direct localization on dtm
        :param row :  line sensor position
        :type row : float
        :param col :  column sensor position
        :type col : float
        :param dtm : dtm model
        :type dtm  : shareloc.dtm
        :return ground position (lon,lat,h)
        :rtype numpy.array
        """
        print("direct localization not yet impelemented for RPC model")
        return None



    def direct_loc_h(self,row,col, alt):
        """evalue loc directe par application du RPC direct"""

        if self.Num_X:
            Xnorm = (col - self.offset_COL)/self.scale_COL
            Ynorm = (row - self.offset_LIG)/self.scale_LIG
            Znorm = (alt - self.offset_ALT)/self.scale_ALT

            if abs(Xnorm)> self.lim_extrapol :
                print("!!!!! l'evaluation au point est extrapolee en colonne ",Xnorm,col)
            if abs(Ynorm)> self.lim_extrapol :
                print("!!!!! l'evaluation au point est extrapolee en ligne ",Ynorm,row)
            if abs(Znorm)> self.lim_extrapol :
                print("!!!!! l'evaluation au point est extrapolee en altitude ",Znorm,alt)

            monomes = array([self.Monomes[i][0]*Xnorm**int(self.Monomes[i][1])*\
                 Ynorm**int(self.Monomes[i][2])*\
                 Znorm**int(self.Monomes[i][3]) for i in range(self.Monomes.__len__())])

            Xout = dot(array(self.Num_X),monomes)/dot(array(self.Den_X),monomes)*self.scale_X+self.offset_X
            Yout = dot(array(self.Num_Y),monomes)/dot(array(self.Den_Y),monomes)*self.scale_Y+self.offset_Y
        else:
            print("les coefficient directs n'ont pas ete definis")
            (Xout,Yout, alt) = (None,None,None)
        return (Xout,Yout, alt)


    def direct_loc_grid_h(self, row0, col0, steprow, stepcol, nbrow, nbcol, alt):
        """calcule une grille de loc directe a partir des RPC directs
         direct localization  grid at constant altitude
         :param row0 :  grid origin (row)
         :type row0 : int
         :param col0 :  grid origin (col)
         :type col0 : int
         :param steprow :  grid step (row)
         :type steprow : int
         :param stepcol :  grid step (col)
         :type stepcol : int
         :param nbrow :  grid nb row
         :type nbrow : int
         :param nbcol :  grid nb col
         :type nbcol : int
         :param alt : altitude of the grid
         :type alt  : float
         :return direct localization grid
         :rtype numpy.array
        """
        gri_lon = zeros((nbrow,nbcol))
        gri_lat = zeros((nbrow,nbcol))
        for c in range(int(nbcol)):
            col = col0 + stepcol*c
            for l in range(int(nbrow)):
                row = row0 + steprow*l
                (gri_lon[l,c],gri_lat[l,c],__) = self.direct_loc_h(row,col,alt)
        return (gri_lon,gri_lat)


    def inverse_loc(self,lon,lat, alt):
        """evalue loc inverse par application du RPC direct"""
        if self.Num_COL:
            Xnorm = (lon - self.offset_X)/self.scale_X
            Ynorm = (lat - self.offset_Y)/self.scale_Y
            Znorm = (alt - self.offset_ALT)/self.scale_ALT

            if abs(Xnorm)> self.lim_extrapol :
                print("!!!!! l'evaluation au point est extrapolee en longitude ",Xnorm,lon)
            if abs(Ynorm)> self.lim_extrapol :
                print("!!!!! l'evaluation au point est extrapolee en latitude ",Ynorm,lat)
            if abs(Znorm)> self.lim_extrapol :
                print("!!!!! l'evaluation au point est extrapolee en altitude ",Znorm,alt)

            monomes = array([self.Monomes[i][0]*Xnorm**int(self.Monomes[i][1])*\
                Ynorm**int(self.Monomes[i][2])*\
                Znorm**int(self.Monomes[i][3]) for i in range(self.Monomes.__len__())])

            Cout = dot(array(self.Num_COL),monomes)/dot(array(self.Den_COL),monomes)*self.scale_COL+self.offset_COL
            Lout = dot(array(self.Num_LIG),monomes)/dot(array(self.Den_LIG),monomes)*self.scale_LIG+self.offset_LIG
        else:
            print("!!!!! les coefficient inverses n'ont pas ete definis")
            (Cout,Lout) = (None,None)
        return (Lout,Cout, True)

    def direct_loc_inverse_iterative(self,row,col,alt,nb_iter_max=10):
        """evalue loc inverse par inversion du RPC inverse        """
        if self.Num_COL:
            #calcul d'une sol approchee: en prend le milieu de la scene
            X = self.offset_X
            Y = self.offset_Y
            (l0,c0, __) = self.inverse_loc(X,Y,alt)

            #precision en pixels
            eps = 1e-6

            k=0
            dc = col - c0
            dl = row - l0
            while abs(dc)>eps and abs(dl)>eps and k<nb_iter_max:
                #evaluer deriv partielles
                (Cdx,Cdy,Ldx,Ldy) = self.calcule_derivees_inv(X,Y,alt)
                det = Cdx*Ldy-Ldx*Cdy
                dX = ( Ldy*dc - Cdy*dl)/det
                dY = (-Ldx*dc + Cdx*dl)/det
                X += dX
                Y += dY
                (l,c, __) = self.inverse_loc(X,Y,alt)
                dc = col - c
                dl = row - l
                k+=1
        else:
            print("!!!!! les coefficient inverses n'ont pas ete definis")
            (X,Y) = (None,None)
        return(X,Y)