# -*- coding: utf-8 -*-
"""
Created on Tue July 28 16:44:35 2020

@author: gresloud
"""
import numpy as np

"""
This module contains the mathematical functions for shareloc
"""

#------------------------------------------------------------------------------
def interpol_bilin(mats,nl,nc,dl,dc):
        """interpole bilineairement une matrice de taille (: , nl,nc)
        au point l,c donnees en indices decimaux de mat
        mats est une liste de mat a interpoler"""
        if (dl < 0):
            i1 = 0
        elif (dl >= nl-1):
            i1 = nl - 2
        else:
            i1 = int(np.floor(dl))
        i2 = i1+1

        if (dc < 0):
            j1 = 0
        elif (dc >= nc-1):
            j1 = nc - 2
        else:
            j1 = int(np.floor(dc))
        j2 = j1+1
        #Coefficients d'interpolation bilineaire

        u = dc - j1
        v = dl - i1
        #Altitude
        matis=[]
        for mat in mats:
            mati = (1-u)*(1-v)*mat[:,i1,j1] + u*(1-v)*mat[:,i1,j2] +\
                   (1-u)*v*mat[:,i2,j1]     + u*v*mat[:,i2,j2]
            matis.append(mati)
        return matis

