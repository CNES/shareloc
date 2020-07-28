# -*- coding: utf-8 -*-
"""
Created on Tue July 28 18:44:35 2020

@author: gresloud
"""
import numpy as np
import os.path as op
import os

#------------------------------------------------------------------------------
def lit_header_hdbabel(fic_mnt_bsq):
    path_header = op.dirname(fic_mnt_bsq)
    name_header = op.basename(fic_mnt_bsq).split('.')[0]+'.hd_babel'
    entete = {}
    hd_babel = {}
    with open(op.join(path_header,name_header),'r') as f:
        txt_header = f.readlines()
    for lig_header in txt_header:
        if lig_header.startswith('>>'):
            ligsplit = lig_header.split('\t')
            entete[ligsplit[1]] = ligsplit[2]

    hd_babel['x0'] = float(entete['LON_REF'])
    hd_babel['y0'] = float(entete['LAT_REF'])
    hd_babel['nc'] = int(entete['NB_LON'])
    hd_babel['nl'] = int(entete['NB_LAT'])
    hd_babel['px'] = float(entete['PAS_LON'])
    hd_babel['py'] = float(entete['PAS_LAT'])
    if int(entete['TYPE_CODE'])==2:
        hd_babel['codage'] = 'int16'
    return hd_babel

#------------------------------------------------------------------------------
def lit_grille_bsq(fic_bsq,nl,nc,codage):
    grille = np.fromfile(fic_bsq,dtype=codage).reshape((nl,nc))
    return grille
#------------------------------------------------------------------------------
def lit_hd_bsq(fic_hd, balises):
    """balises est un dico
    dico_out['nom'] = ('nom_bsq',type)
    """
    dico_out = {}
    with open(op.join(fic_hd),'r') as f:
        txt_header = f.readlines()
    entete = {}

    for i in range(0,txt_header.__len__(),2):
        entete[txt_header[i].strip()] = txt_header[i+1].strip()

    for var in balises:
        (nom, form) = balises[var]
        dico_out[var] = form(entete[nom])

    return dico_out
