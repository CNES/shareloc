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


class FonctRatD:
	def __init__(self,fichier_dimap_or_euclide_d,fichier_euclide_i=''):

		self.offset_COL	= None
		self.scale_COL	= None
		self.offset_LIG	= None
		self.scale_LIG	= None
		self.offset_ALT	= None
		self.scale_ALT	= None
		self.offset_X	= None
		self.scale_X	= None
		self.offset_Y	= None
		self.scale_Y	= None
		self.Monomes	= None
		self.Num_X		= None
		self.Den_X		= None
		self.Num_Y		= None
		self.Den_Y		= None
		self.Num_COL	= None
		self.Den_COL	= None
		self.Num_LIG	= None
		self.Den_LIG	= None

		ordre_monomes_LAI = \
				[[0,0,0],[1,0,0],[0,1,0],\
				[0,0,1],[1,1,0],[1,0,1],\
				[0,1,1],[2,0,0],[0,2,0],\
				[0,0,2],[1,1,1],[3,0,0],\
				[1,2,0],[1,0,2],[2,1,0],\
				[0,3,0],[0,1,2],[2,0,1],\
				[0,2,1],[0,0,3]]

		self.Monomes	= ordre_monomes_LAI

		if basename(fichier_dimap_or_euclide_d).endswith('XML'):
			xmldoc= minidom.parse(fichier_dimap_or_euclide_d)
			GLOBAL_RFM	= xmldoc.getElementsByTagName('Global_RFM')
			RFM_Validity 	= xmldoc.getElementsByTagName('RFM_Validity')
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


			self.offset_COL	= B_col
			self.scale_COL	= A_col
			self.offset_LIG	= B_row
			self.scale_LIG	= A_row
			self.offset_ALT	= B_alt
			self.scale_ALT	= A_alt
			self.offset_X	= B_lon
			self.scale_X	= A_lon
			self.offset_Y	= B_lat
			self.scale_Y	= A_lat

			self.Num_X		= coeff_LON[0:20]
			self.Den_X		= coeff_LON[20::]
			self.Num_Y		= coeff_LAT[0:20]
			self.Den_Y		= coeff_LAT[20::]
			self.Num_COL	= coeff_COL[0:20]
			self.Den_COL	= coeff_COL[20::]
			self.Num_LIG	= coeff_LIG[0:20]
			self.Den_LIG	= coeff_LIG[20::]



		else:
			fid = open(fichier_dimap_or_euclide_d,'r')
			txt = fid.readlines()
			fid.close()
			lsep = renvoie_linesep(txt)

			ind_debut_PX = txt.index('>>\tCOEFF POLYNOME PXOUT'+lsep)
			ind_debut_QX = txt.index('>>\tCOEFF POLYNOME QXOUT'+lsep)
			ind_debut_PY = txt.index('>>\tCOEFF POLYNOME PYOUT'+lsep)
			ind_debut_QY = txt.index('>>\tCOEFF POLYNOME QYOUT'+lsep)


			coeff_PX_str = txt[ind_debut_PX+1:ind_debut_PX+21]
			coeff_QX_str = txt[ind_debut_QX+1:ind_debut_QX+21]
			coeff_PY_str = txt[ind_debut_PY+1:ind_debut_PY+21]
			coeff_QY_str = txt[ind_debut_QY+1:ind_debut_QY+21]

			coeff_PX_d = [float(coeff.split()[1]) for coeff in coeff_PX_str]
			coeff_QX_d = [float(coeff.split()[1]) for coeff in coeff_QX_str]
			coeff_PY_d = [float(coeff.split()[1]) for coeff in coeff_PY_str]
			coeff_QY_d = [float(coeff.split()[1]) for coeff in coeff_QY_str]

			for l in txt:
				if l.startswith('>>\tTYPE_OBJET'):
					if l.split()[-1].endswith('Inverse'):
						type_fic = 'I'
					if l.split()[-1].endswith('Directe'):
						type_fic = 'D'
				if l.startswith('>>\tXIN_OFFSET'):
					lsplit = l.split()
					(offset_COL_d,scale_COL_d) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tYIN_OFFSET'):
					lsplit = l.split()
					(offset_LIG_d,scale_LIG_d) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tZIN_OFFSET'):
					lsplit = l.split()
					(offset_ALT_d,scale_ALT_d) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tXOUT_OFFSET'):
					lsplit = l.split()
					(offset_X_d,scale_X_d) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tYOUT_OFFSET'):
					lsplit = l.split()
					(offset_Y_d,scale_Y_d) = (float(lsplit[4]),float(lsplit[5]))
			if type_fic != 'D':
				print("le fichier euclide direct est inconherent (le sens est marque Inverse)")


			message = "!!!!! les fichiers directs et inverses sont incoherents car {} differents"
			if not self.offset_COL:
				self.offset_COL = offset_COL_d
			elif self.offset_COL != offset_COL_d:
				print(message.format("offset_COL"))

			if not self.scale_COL:
				self.scale_COL = scale_COL_d
			elif self.scale_COL != scale_COL_d:
				print(message.format("scale_COL"))

			if not self.offset_LIG:
				self.offset_LIG = offset_LIG_d
			elif self.offset_LIG != offset_LIG_d:
				print(message.format("offset_LIG"))

			if not self.scale_LIG:
				self.scale_LIG = scale_LIG_d
			elif self.scale_LIG != scale_LIG_d:
				pprint(message.format("scale_LIG"))

			if not self.offset_ALT:
				self.offset_ALT = offset_ALT_d
			elif self.offset_ALT != offset_ALT_d:
				print(message.format("offset_ALT"))

			if not self.scale_ALT:
				self.scale_ALT = scale_ALT_d
			elif self.scale_ALT != scale_ALT_d:
				print(message.format("scale_ALT"))

			if not self.offset_X:
				self.offset_X = offset_X_d
			elif self.offset_X != offset_X_d:
				print(message.format("offset_X"))

			if not self.scale_X:
				self.scale_X = scale_X_d
			elif self.scale_X != scale_X_d:
				print(message.format("scale_X"))

			if not self.offset_Y:
				self.offset_Y = offset_Y_d
			elif self.offset_Y != offset_Y_d:
				print(message.format("offset_Y"))

			if not self.scale_Y:
				self.scale_Y = scale_Y_d
			elif self.scale_Y != scale_Y_d:
				print(message.format("scale_Y"))

			self.Num_X	= coeff_PX_d
			self.Den_X	= coeff_QX_d
			self.Num_Y	= coeff_PY_d
			self.Den_Y	= coeff_QY_d

		if fichier_euclide_i:
			fid = open(fichier_euclide_i,'r')
			txt = fid.readlines()
			fid.close()
			lsep = renvoie_linesep(txt)

			ind_debut_PX = txt.index('>>\tCOEFF POLYNOME PXOUT'+lsep)
			ind_debut_QX = txt.index('>>\tCOEFF POLYNOME QXOUT'+lsep)
			ind_debut_PY = txt.index('>>\tCOEFF POLYNOME PYOUT'+lsep)
			ind_debut_QY = txt.index('>>\tCOEFF POLYNOME QYOUT'+lsep)

			coeff_PX_str = txt[ind_debut_PX+1:ind_debut_PX+21]
			coeff_QX_str = txt[ind_debut_QX+1:ind_debut_QX+21]
			coeff_PY_str = txt[ind_debut_PY+1:ind_debut_PY+21]
			coeff_QY_str = txt[ind_debut_QY+1:ind_debut_QY+21]

			coeff_PX_i = [float(coeff.split()[1]) for coeff in coeff_PX_str]
			coeff_QX_i = [float(coeff.split()[1]) for coeff in coeff_QX_str]
			coeff_PY_i = [float(coeff.split()[1]) for coeff in coeff_PY_str]
			coeff_QY_i = [float(coeff.split()[1]) for coeff in coeff_QY_str]

			for l in txt:
				if l.startswith('>>\tTYPE_OBJET'):
					if l.split()[-1].endswith('Inverse'):
						type_fic = 'I'
					if l.split()[-1].endswith('Directe'):
						type_fic = 'D'
				if l.startswith('>>\tXIN_OFFSET'):
					lsplit = l.split()
					(offset_X_i,scale_X_i) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tYIN_OFFSET'):
					lsplit = l.split()
					(offset_Y_i,scale_Y_i) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tZIN_OFFSET'):
					lsplit = l.split()
					(offset_ALT_i,scale_ALT_i) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tXOUT_OFFSET'):
					lsplit = l.split()
					(offset_COL_i,scale_COL_i) = (float(lsplit[4]),float(lsplit[5]))
				if l.startswith('>>\tYOUT_OFFSET'):
					lsplit = l.split()
					(offset_LIG_i,scale_LIG_i) = (float(lsplit[4]),float(lsplit[5]))

			if type_fic != 'I':
				print("!!!!! le fichier euclide inverse est inconherent (le sens est marque Direct)")


			if not self.offset_COL:
				self.offset_COL = offset_COL_i
			elif self.offset_COL != offset_COL_i:
				print(message.format("offset_COL"))
			if not self.scale_COL:
				self.scale_COL = scale_COL_i
			elif self.scale_COL != scale_COL_i:
				print(message.format("scale_COL"))

			if not self.offset_LIG:
				self.offset_LIG = offset_LIG_i
			elif self.offset_LIG != offset_LIG_i:
				print(message.format("offset_LIG"))

			if not self.scale_LIG:
				self.scale_LIG = scale_LIG_i
			elif self.scale_LIG != scale_LIG_i:
				print(message.format("scale_LIG"))

			if not self.offset_ALT:
				self.offset_ALT = offset_ALT_i
			elif self.offset_ALT != offset_ALT_i:
				print(message.format("offset_ALT"))

			if not self.scale_ALT:
				self.scale_ALT = scale_ALT_i
			elif self.scale_ALT != scale_ALT_i:
				print(message.format("scale_ALT"))

			if not self.offset_X:
				self.offset_X = offset_X_i
			elif self.offset_X != offset_X_i:
				print(message.format("offset_X"))

			if not self.scale_X:
				self.scale_X = scale_X_i
			elif self.scale_X != scale_X_i:
				print(message.format("scale_X"))

			if not self.offset_Y:
				self.offset_Y = offset_Y_i
			elif self.offset_Y != offset_Y_i:
				print(message.format("offset_Y"))

			if not self.scale_Y:
				self.scale_Y = scale_Y_i
			elif self.scale_Y != scale_Y_i:
				print(message.format("scale_Y"))

			self.Num_COL	= coeff_PX_i
			self.Den_COL	= coeff_QX_i
			self.Num_LIG	= coeff_PY_i
			self.Den_LIG	= coeff_QY_i

	def direct_loc_h(self,col,lig, alt):
		if self.Num_X:
			Xnorm = (col - self.offset_COL)/self.scale_COL
			Ynorm = (lig - self.offset_LIG)/self.scale_LIG
			Znorm = (alt - self.offset_ALT)/self.scale_ALT

			message =  "!!!!! l'evaluation au point est extrapolee en {} {} {}"
			if abs(Xnorm)> 1.001 :
				print(message.format("colonne", Xnorm, col))
			if abs(Ynorm)> 1.001 :
				print(message.format("ligne", Ynorm, lig))
			if abs(Znorm)> 1.001 :
				print(message.format("altitude", Znorm, alt))

			monomes = array([Xnorm**int(self.Monomes[i][0])*\
				Ynorm**int(self.Monomes[i][1])*\
				Znorm**int(self.Monomes[i][2]) for i in range(self.Monomes.__len__())])
			Xout = dot(array(self.Num_X),monomes)/dot(array(self.Den_X),monomes)*self.scale_X+self.offset_X
			Yout = dot(array(self.Num_Y),monomes)/dot(array(self.Den_Y),monomes)*self.scale_Y+self.offset_Y
		else:
			print("les coefficient directs n'ont pas ete definis")
			(Xout,Yout) = (None,None)
		return (Xout,Yout)

	def direct_loc_grid_h(self,c0,l0,pascol,paslig,nbcol,nblig, alt):
		gri_lon = zeros((nblig,nbcol))
		gri_lat = zeros((nblig,nbcol))
		for c in range(int(nbcol)):
			col = c0 + pascol*c
			for l in range(int(nblig)):
				lig = l0 + paslig*l
				(gri_lon[l,c],gri_lat[l,c]) = self.evalue_loc_d(col,lig,alt)
		return (gri_lon,gri_lat)

	def inverse_loc(self,lon,lat, alt):
		if self.Num_COL:
			Xnorm = (lon - self.offset_X)/self.scale_X
			Ynorm = (lat - self.offset_Y)/self.scale_Y
			Znorm = (alt - self.offset_ALT)/self.scale_ALT

			message = "!!!!! l'evaluation au point est extrapolee en {} {} {}"
			if abs(Xnorm)> 1.001 :
				print(message.format("longitude", Xnorm, lon))
			if abs(Ynorm)> 1.001 :
				print(message.format("latitude", Ynorm, lat))
			if abs(Znorm)> 1.001 :
				print(message.format("altitude", Znorm, alt))



			monomes = array([Xnorm**int(self.Monomes[i][0])*\
				Ynorm**int(self.Monomes[i][1])*\
				Znorm**int(self.Monomes[i][2]) for i in range(self.Monomes.__len__())])

			Cout = dot(array(self.Num_COL),monomes)/dot(array(self.Den_COL),monomes)*self.scale_COL+self.offset_COL
			Lout = dot(array(self.Num_LIG),monomes)/dot(array(self.Den_LIG),monomes)*self.scale_LIG+self.offset_LIG
		else:
			print("!!!!! les coefficient inverses n'ont pas ete definis")
			(Cout,Lout) = (None,None)
		return (Cout,Lout)


class FonctRatD_simple:
	def __init__(self,Kc,Kl,Klon,Klat,B_col,A_col,B_row,A_row,B_lon,A_lon,B_lat,A_lat,B_alt,A_alt):

		self.offset_COL	= None
		self.scale_COL	= None
		self.offset_LIG	= None
		self.scale_LIG	= None
		self.offset_ALT	= None
		self.scale_ALT	= None
		self.offset_X	= None
		self.scale_X	= None
		self.offset_Y	= None
		self.scale_Y	= None
		self.Monomes	= None
		self.Num_X		= None
		self.Den_X		= None
		self.Num_Y		= None
		self.Den_Y		= None
		self.Num_COL	= None
		self.Den_COL	= None
		self.Num_LIG	= None
		self.Den_LIG	= None

		ordre_monomes_LAI = \
				[[0,0,0],[1,0,0],[0,1,0],\
				[0,0,1],[1,1,0],[1,0,1],\
				[0,1,1],[2,0,0],[0,2,0],\
				[0,0,2],[1,1,1],[3,0,0],\
				[1,2,0],[1,0,2],[2,1,0],\
				[0,3,0],[0,1,2],[2,0,1],\
				[0,2,1],[0,0,3]]

		self.Monomes	= ordre_monomes_LAI

		self.offset_COL	= B_col
		self.scale_COL	= A_col
		self.offset_LIG	= B_row
		self.scale_LIG	= A_row
		self.offset_ALT	= B_alt
		self.scale_ALT	= A_alt
		self.offset_X	= B_lon
		self.scale_X	= A_lon
		self.offset_Y	= B_lat
		self.scale_Y	= A_lat

		self.Num_COL	= Kc[0:20]
		self.Den_COL	= Kc[20::]
		self.Num_LIG	= Kl[0:20]
		self.Den_LIG	= Kl[20::]
		self.Num_X		= Klon[0:20]
		self.Den_X		= Klon[20::]
		self.Num_Y		= Klat[0:20]
		self.Den_Y		= Klat[20::]


	def direct_loc_h(self,col,lig, alt):
		if self.Num_X:
			Xnorm = (col - self.offset_COL)/self.scale_COL
			Ynorm = (lig - self.offset_LIG)/self.scale_LIG
			Znorm = (alt - self.offset_ALT)/self.scale_ALT

			message = "!!!!! l'evaluation au point est extrapolee en {} {} {}"
			if abs(Xnorm)> 1.001 :
				print(message.format("colonne", Xnorm, col))
			if abs(Ynorm)> 1.001 :
				print(message.format("ligne", Ynorm, lig))
			if abs(Znorm)> 1.001 :
				print(message.format("altitude", Znorm, alt))

			monomes = array([Xnorm**int(self.Monomes[i][0])*\
				Ynorm**int(self.Monomes[i][1])*\
				Znorm**int(self.Monomes[i][2]) for i in range(self.Monomes.__len__())])
			Xout = dot(array(self.Num_X),monomes)/dot(array(self.Den_X),monomes)*self.scale_X+self.offset_X
			Yout = dot(array(self.Num_Y),monomes)/dot(array(self.Den_Y),monomes)*self.scale_Y+self.offset_Y
		else:
			print("les coefficient directs n'ont pas ete definis")
			(Xout,Yout) = (None,None)
		return (Xout,Yout)

	def direct_loc_grid_h(self,c0,l0,pascol,paslig,nbcol,nblig, alt):
		gri_lon = zeros((nblig,nbcol))
		gri_lat = zeros((nblig,nbcol))
		for c in range(int(nbcol)):
			col = c0 + pascol*c
			for l in range(int(nblig)):
				lig = l0 + paslig*l
				(gri_lon[l,c],gri_lat[l,c]) = self.evalue_loc_d(col,lig,alt)
		return (gri_lon,gri_lat)

	def inverse_loc(self,lon,lat, alt):
		if self.Num_COL != None:
			Xnorm = (lon - self.offset_X)/self.scale_X
			Ynorm = (lat - self.offset_Y)/self.scale_Y
			Znorm = (alt - self.offset_ALT)/self.scale_ALT

			message = "!!!!! l'evaluation au point est extrapolee en {} {} {}"
			if abs(Xnorm)> 1.001 :
				print(message.format("longitude", Xnorm, lon))
			if abs(Ynorm)> 1.001 :
				print(message.format("latitude", Ynorm, lat))
			if abs(Znorm)> 1.001 :
				print(message.format("altitude", Znorm, alt))


			monomes = array([Xnorm**int(self.Monomes[i][0])*\
				Ynorm**int(self.Monomes[i][1])*\
				Znorm**int(self.Monomes[i][2]) for i in range(self.Monomes.__len__())])

			Cout = dot(array(self.Num_COL),monomes)/dot(array(self.Den_COL),monomes)*self.scale_COL+self.offset_COL
			Lout = dot(array(self.Num_LIG),monomes)/dot(array(self.Den_LIG),monomes)*self.scale_LIG+self.offset_LIG
		else:
			print("!!!!! les coefficient inverses n'ont pas ete definis")
			(Cout,Lout) = (None,None)
		return (Cout,Lout)
