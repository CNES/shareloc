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

import rasterio as rio
from xml.dom import minidom
from os.path import basename
import numpy as np
from numba import njit, prange


def renvoie_linesep(txt_liste_lines):
	"""Renvoie le separateur de ligne d'un texte sous forme de liste de lignes
	Obtenu par readlines
	"""
	if txt_liste_lines[0].endswith('\r\n'):
		line_sep = '\r\n'
	elif txt_liste_lines[0].endswith('\n'):
		line_sep = '\n'
	return line_sep


def parse_coeff_line(coeff_str):
    """
    split str coef to float list
    :param coeff_str : line coef
    :type coeff_str : str
    :return coeff list
    :rtype list()
    """
    return [float(el) for el in coeff_str.split()]

def identify_dimap(xml_file):
    """
    parse xml file to identify dimap and its version
    :param xml_file : dimap rpc file
    :type xml_file : str
    :return dimap info : dimap_version and None if not an dimap file
    :rtype str
    """
    try :
        xmldoc = minidom.parse(xml_file)
        mtd = xmldoc.getElementsByTagName('Metadata_Identification')
        mtd_format = mtd[0].getElementsByTagName('METADATA_FORMAT')[0].firstChild.data
        if mtd_format == 'DIMAP_PHR':
            version_tag = 'METADATA_PROFILE'
        else:
            version_tag = 'METADATA_FORMAT'
        version = mtd[0].getElementsByTagName(version_tag)[0].attributes.items()[0][1]
        return version
    except:
        return None


def identify_euclidium_rpc(eucl_file):
    """
    parse file to identify if it is an euclidium model (starts with '>>')
    :param eucl_file : euclidium rpc file
    :type eucl_file : str
    :return  True if euclidium rpc has been identified, False otherwise
    :rtype Boolean
    """
    try :
        with open(eucl_file) as f:
                content = f.readlines()
        is_eucl = True
        for line in content:
            if not line.startswith('>>') and line !='\n':
                is_eucl = False
            return is_eucl
    except:
        return False

def identify_ossim_kwl(ossim_kwl_file):
    """
    parse geom file to identify if it is an ossim model
    :param ossim_kwl_file : ossim keyword list file
    :type ossim_kwl_file : str
    :return ossim kwl info : ossimmodel or None if not an ossim kwl file
    :rtype str
    """
    try :
        with open(ossim_kwl_file) as f:
            content = f.readlines()

        geom_dict = dict()
        for line in content:
            (key, val) = line.split(': ')
            geom_dict[key] = val.rstrip()
        if 'type' in geom_dict.keys():
            if geom_dict['type'].strip().startswith('ossim') :
                return geom_dict['type'].strip()
            else:
                return None
    except:
        return None

def identify_geotiff_rpc(image_filename):
    """
    read image file to identify if it is a geotiff which contains RPCs
    :param image_filename : image_filename
    :type image_filename : str
    :return rpc info : rpc dict or None  if not a geotiff with rpc
    :rtype str
    """
    try :
        dataset = rio.open(image_filename)
        rpc_dict = dataset.tags(ns='RPC')
        if not rpc_dict:
            return None
        else:
            return rpc_dict
    except:
        return None


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

    for line in txt:
        if line.startswith('>>\tTYPE_OBJET'):
            if line.split()[-1].endswith('Inverse'):
                parsed_file['type_fic'] = 'I'
            if line.split()[-1].endswith('Directe'):
                parsed_file['type_fic'] = 'D'

    lsep = renvoie_linesep(txt)

    ind_debut_PX = txt.index('>>\tCOEFF POLYNOME PXOUT' + lsep)
    ind_debut_QX = txt.index('>>\tCOEFF POLYNOME QXOUT' + lsep)
    ind_debut_PY = txt.index('>>\tCOEFF POLYNOME PYOUT' + lsep)
    ind_debut_QY = txt.index('>>\tCOEFF POLYNOME QYOUT' + lsep)

    coeff_PX_str = txt[ind_debut_PX + 1:ind_debut_PX + 21]
    coeff_QX_str = txt[ind_debut_QX + 1:ind_debut_QX + 21]
    coeff_PY_str = txt[ind_debut_PY + 1:ind_debut_PY + 21]
    coeff_QY_str = txt[ind_debut_QY + 1:ind_debut_QY + 21]

    poly_coeffs = dict()
    if parsed_file['type_fic'] == 'I':
        poly_coeffs['Num_COL'] = [float(coeff.split()[1]) for coeff in coeff_PX_str]
        poly_coeffs['Den_COL'] = [float(coeff.split()[1]) for coeff in coeff_QX_str]
        poly_coeffs['Num_LIG'] = [float(coeff.split()[1]) for coeff in coeff_PY_str]
        poly_coeffs['Den_LIG'] = [float(coeff.split()[1]) for coeff in coeff_QY_str]
    else:
        poly_coeffs['Num_X'] = [float(coeff.split()[1]) for coeff in coeff_PX_str]
        poly_coeffs['Den_X'] = [float(coeff.split()[1]) for coeff in coeff_QX_str]
        poly_coeffs['Num_Y'] = [float(coeff.split()[1]) for coeff in coeff_PY_str]
        poly_coeffs['Den_Y'] = [float(coeff.split()[1]) for coeff in coeff_QY_str]

    parsed_file['poly_coeffs'] = poly_coeffs
    # list [offset , scale]
    normalisation_coeff = dict()
    for l in txt:
        if l.startswith('>>\tXIN_OFFSET'):
            lsplit = l.split()
            if parsed_file['type_fic'] == 'I':
                param ='X'
            else:
                param = 'COL'
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tYIN_OFFSET'):
            if parsed_file['type_fic'] == 'I':
                param ='Y'
            else:
                param = 'LIG'
            lsplit = l.split()
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tZIN_OFFSET'):
            lsplit = l.split()
            normalisation_coeff['ALT'] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tXOUT_OFFSET'):
            lsplit = l.split()
            if parsed_file['type_fic'] == 'D':
                param ='X'
            else:
                param = 'COL'
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
        if l.startswith('>>\tYOUT_OFFSET'):
            lsplit = l.split()
            if parsed_file['type_fic'] == 'D':
                param ='Y'
            else:
                param = 'LIG'
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
    parsed_file['normalisation_coeffs'] = normalisation_coeff
    return parsed_file


def check_coeff_consistency(dict1, dict2):
    """
    print an error message inf normalisations coeff are not consistent
    :param dict1 : normalisation coeffs 1
    :type dict1 : dict
    :param dict2 : normalisation coeffs 2
    :type dict2 : dict

    """
    for key, value in dict1.items() :
        if dict2[key] != value:
            print("normalisation coeffs are different between"
                  " direct en inverse one : {} : {} {}".format(key,value,dict2[key]))

class RPC:
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

        self.Monomes    = np.array(ordre_monomes_LAI)

        #coefficient des degres monomes avec derivation 1ere variable
        self.monomes_deriv_1 = np.array(
                [[0,0,0,0],[1,0,0,0],[0,0,1,0],
                 [0,0,0,1],[1,0,1,0],[1,0,0,1],
                 [0,0,1,1],[2,1,0,0],[0,0,2,0],
                 [0,0,0,2],[1,0,1,1],[3,2,0,0],
                 [1,0,2,0],[1,0,0,2],[2,1,1,0],
                 [0,0,3,0],[0,0,1,2],[2,1,0,1],
                 [0,0,2,1],[0,0,0,3]])

        #coefficient des degres monomes avec derivation 1ere variable
        self.monomes_deriv_2 =  np.array(
                [[0,0,0,0],[0,1,0,0],[1,0,0,0],
                 [0,0,0,1],[1,1,0,0],[0,1,0,1],
                 [1,0,0,1],[0,2,0,0],[2,0,1,0],
                 [0,0,0,2],[1,1,0,1],[0,3,0,0],
                 [2,1,1,0],[0,1,0,2],[1,2,0,0],
                 [3,0,2,0],[1,0,0,2],[0,2,0,1],
                 [2,0,1,1],[0,0,0,3]])

    @classmethod
    def from_dimap(cls, dimap_filepath, topleftconvention=True):
        """ load from dimap
        param dimap_filepath  : dimap xml file
    	:type dimap_filepath  : str
        :param topleftconvention  : [0,0] position
    	:type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        dimap_version = identify_dimap(dimap_filepath)
        if identify_dimap(dimap_filepath) is not None:
            if float(dimap_version) < 2.0:
                return cls.from_dimap_v1(dimap_filepath, topleftconvention)
            else:
                return cls.from_dimap_v2(dimap_filepath, topleftconvention)
        else:
            ValueError("can''t read dimap file")
            return None

    @classmethod
    def from_dimap_v2(cls, dimap_filepath, topleftconvention=True):
        """ load from dimap  v2
        :param dimap_filepath  : dimap xml file
    	:type dimap_filepath  : str
        :param topleftconvention  : [0,0] position
    	:type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """

        rpc_params = dict()

        if not basename(dimap_filepath).endswith('XML'.upper()):
            raise ValueError("dimap must ends with .xml")

        xmldoc = minidom.parse(dimap_filepath)

        mtd = xmldoc.getElementsByTagName('Metadata_Identification')
        version = mtd[0].getElementsByTagName('METADATA_FORMAT')[0].attributes.items()[0][1]
        rpc_params['driver_type'] = 'dimap_v' + version
        GLOBAL_RFM    = xmldoc.getElementsByTagName('Global_RFM')[0]
        direct_coeffs = GLOBAL_RFM.getElementsByTagName('Direct_Model')[0]
        rpc_params['Num_X'] = [float(direct_coeffs.getElementsByTagName('SAMP_NUM_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        rpc_params['Den_X'] = [float(direct_coeffs.getElementsByTagName('SAMP_DEN_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        rpc_params['Num_Y'] = [float(direct_coeffs.getElementsByTagName('LINE_NUM_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        rpc_params['Den_Y'] = [float(direct_coeffs.getElementsByTagName('LINE_DEN_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        inverse_coeffs = GLOBAL_RFM.getElementsByTagName('Inverse_Model')[0]
        rpc_params['Num_COL'] = [float(inverse_coeffs.getElementsByTagName('SAMP_NUM_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        rpc_params['Den_COL'] = [float(inverse_coeffs.getElementsByTagName('SAMP_DEN_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        rpc_params['Num_LIG'] = [float(inverse_coeffs.getElementsByTagName('LINE_NUM_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        rpc_params['Den_LIG'] = [float(inverse_coeffs.getElementsByTagName('LINE_DEN_COEFF_{}'.
                                                                          format(index))[0].firstChild.data)
                for index in range(1,21)]
        normalisation_coeffs = GLOBAL_RFM.getElementsByTagName('RFM_Validity')[0]
        rpc_params['offset_COL']    = float(normalisation_coeffs.getElementsByTagName('SAMP_OFF')[0].firstChild.data)
        rpc_params['scale_COL']    = float(normalisation_coeffs.getElementsByTagName('SAMP_SCALE')[0].firstChild.data)
        rpc_params['offset_LIG']    = float(normalisation_coeffs.getElementsByTagName('LINE_OFF')[0].firstChild.data)
        rpc_params['scale_LIG']    = float(normalisation_coeffs.getElementsByTagName('LINE_SCALE')[0].firstChild.data)
        rpc_params['offset_ALT']    = float(normalisation_coeffs.getElementsByTagName('HEIGHT_OFF')[0].firstChild.data)
        rpc_params['scale_ALT']    = float(normalisation_coeffs.getElementsByTagName('HEIGHT_SCALE')[0].firstChild.data)
        rpc_params['offset_X']    = float(normalisation_coeffs.getElementsByTagName('LONG_OFF')[0].firstChild.data)
        rpc_params['scale_X']    = float(normalisation_coeffs.getElementsByTagName('LONG_SCALE')[0].firstChild.data)
        rpc_params['offset_Y']    = float(normalisation_coeffs.getElementsByTagName('LAT_OFF')[0].firstChild.data)
        rpc_params['scale_Y']    = float(normalisation_coeffs.getElementsByTagName('LAT_SCALE')[0].firstChild.data)
        rpc_params['offset_COL'] -= 1.5
        rpc_params['offset_LIG'] -= 1.5
        #If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params['offset_COL'] += 0.5
            rpc_params['offset_LIG'] += 0.5
        return cls(rpc_params)



    @classmethod
    def from_dimap_v1(cls, dimap_filepath, topleftconvention=True):
        """ load from dimap  v1
        :param dimap_filepath  : dimap xml file
	    :type dimap_filepath  : str
        :param topleftconvention  : [0,0] position
	    :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel 
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        rpc_params = dict()

        if not basename(dimap_filepath).endswith('XML'.upper()):
            raise ValueError("dimap must ends with .xml")


        xmldoc= minidom.parse(dimap_filepath)

        mtd = xmldoc.getElementsByTagName('Metadata_Identification')
        version = mtd[0].getElementsByTagName('METADATA_PROFILE')[0].attributes.items()[0][1]
        rpc_params['driver_type'] = 'dimap_v' + version

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
        rpc_params['offset_COL'] -= 1.5
        rpc_params['offset_LIG'] -= 1.5
        #If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params['offset_COL'] += 0.5
            rpc_params['offset_LIG'] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_geotiff(cls, image_filename, topleftconvention=True):
        """ Load from a  geotiff image file
        :param image_filename  : image filename
    	:type image_filename  : str
        :param topleftconvention  : [0,0] position
    	:type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        dataset = rio.open(image_filename)
        rpc_dict = dataset.tags(ns='RPC')
        if not rpc_dict:
            print("{} doesn't contains RPCS ".format(image_filename))
            raise ValueError
        rpc_params = dict()
        rpc_params['Den_LIG'] = parse_coeff_line(rpc_dict['LINE_DEN_COEFF'])
        rpc_params['Num_LIG'] = parse_coeff_line(rpc_dict['LINE_NUM_COEFF'])
        rpc_params['Num_COL'] = parse_coeff_line(rpc_dict['SAMP_NUM_COEFF'])
        rpc_params['Den_COL'] = parse_coeff_line(rpc_dict['SAMP_DEN_COEFF'])
        rpc_params['offset_COL']   = float(rpc_dict["SAMP_OFF"])
        rpc_params['scale_COL']    = float(rpc_dict["SAMP_SCALE"])
        rpc_params['offset_LIG']   = float(rpc_dict["LINE_OFF"])
        rpc_params['scale_LIG']    = float(rpc_dict["LINE_SCALE"])
        rpc_params['offset_ALT']   = float(rpc_dict["HEIGHT_OFF"])
        rpc_params['scale_ALT']    = float(rpc_dict["HEIGHT_SCALE"])
        rpc_params['offset_X']   = float(rpc_dict["LONG_OFF"])
        rpc_params['scale_X']    = float(rpc_dict["LONG_SCALE"])
        rpc_params['offset_Y']   = float(rpc_dict["LAT_OFF"])
        rpc_params['scale_Y']    = float(rpc_dict["LAT_SCALE"])
        #inverse coeff are not defined
        rpc_params['Num_X'] = None
        rpc_params['Den_X'] = None
        rpc_params['Num_Y'] = None
        rpc_params['Den_Y'] = None
        #If top left convention, 0.5 pixel shift added on col/row offsets
        rpc_params['offset_COL'] -= 0.5
        rpc_params['offset_LIG'] -= 0.5
        if topleftconvention:
            rpc_params['offset_COL'] += 0.5
            rpc_params['offset_LIG'] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_ossim_kwl(cls, ossim_kwl_filename, topleftconvention=True):
        """ Load from a geom file
        :param topleftconvention  : [0,0] position
	    :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel 
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        rpc_params = dict()
        #OSSIM keyword list
        rpc_params['driver_type'] = 'ossim_kwl'
        with open(ossim_kwl_filename) as f:
            content = f.readlines()

        geom_dict = dict()
        for line in content:
            (key, val) = line.split(': ')
            geom_dict[key] = val.rstrip()

        rpc_params['Den_LIG']= [np.nan] * 20
        rpc_params['Num_LIG'] = [np.nan] * 20
        rpc_params['Den_COL']= [np.nan] * 20
        rpc_params['Num_COL'] = [np.nan] * 20
        for index in range(0, 20):
            axis = "line"
            num_den = "den"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params['Den_LIG'][index] = float(geom_dict[key])
            num_den = "num"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params['Num_LIG'][index] = float(geom_dict[key])
            axis = "samp"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params['Num_COL'][index] = float(geom_dict[key])
            num_den = "den"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params['Den_COL'][index] = float(geom_dict[key])
        rpc_params['offset_COL']    = float(geom_dict["samp_off"])
        rpc_params['scale_COL']    = float(geom_dict["samp_scale"])
        rpc_params['offset_LIG']    = float(geom_dict["line_off"])
        rpc_params['scale_LIG']    = float(geom_dict["line_scale"])
        rpc_params['offset_ALT']    = float(geom_dict["height_off"])
        rpc_params['scale_ALT']    = float(geom_dict["height_scale"])
        rpc_params['offset_X']    = float(geom_dict["long_off"])
        rpc_params['scale_X']    = float(geom_dict["long_scale"])
        rpc_params['offset_Y']    = float(geom_dict["lat_off"])
        rpc_params['scale_Y']    = float(geom_dict["lat_scale"])
        #inverse coeff are not defined
        rpc_params['Num_X'] = None
        rpc_params['Den_X'] = None
        rpc_params['Num_Y'] = None
        rpc_params['Den_Y'] = None
        #If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params['offset_COL'] += 0.5
            rpc_params['offset_LIG'] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_euclidium(cls, primary_euclidium_coeff, secondary_euclidium_coeff=None, topleftconvention=True):
        """ load from euclidium
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        :param primary_euclidium_coeff  : primary euclidium coefficients file (can be either direct or inverse)
        :type primary_euclidium_coeff  : str
        :param secondary_euclidium_coeff  : optional secondary euclidium coeff coefficients file
            (can be either direct or inverse)
        :type secondary_euclidium_coeff  : str
        If False : [0,0] is at the center of the Top Left pixel 
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        rpc_params = dict()
        rpc_params['driver_type'] = 'euclidium'

        #lecture fichier euclidium
        primary_coeffs = read_eucl_file(primary_euclidium_coeff)

        # info log
        print("primary euclidium file is of {} type".format(primary_coeffs['type_fic']))

        rpc_params['Num_X'] = None
        rpc_params['Den_X'] = None
        rpc_params['Num_Y'] = None
        rpc_params['Den_Y'] = None
        rpc_params['Num_COL'] = None
        rpc_params['Den_COL'] = None
        rpc_params['Num_LIG'] = None
        rpc_params['Den_LIG'] = None

        for key, value in primary_coeffs['poly_coeffs'].items():
            rpc_params[key] = value

        rpc_params['normalisation_coeffs'] = primary_coeffs['normalisation_coeffs']
        for key, value in primary_coeffs['normalisation_coeffs'].items():
            rpc_params['offset_' + key] = value[0]
            rpc_params['scale_' + key] = value[1]

        if secondary_euclidium_coeff is not None :
            secondary_coeffs = read_eucl_file(secondary_euclidium_coeff)
            print("secondary euclidium file is of {} type".format(secondary_coeffs['type_fic']))
            check_coeff_consistency(primary_coeffs['normalisation_coeffs'], secondary_coeffs['normalisation_coeffs'])

            for key, value in secondary_coeffs['poly_coeffs'].items():
                rpc_params[key] = value

        rpc_params['offset_COL'] -= 1
        rpc_params['offset_LIG'] -= 1
        #If top left convention, 0.5 pixel shift added on col/row offsets

        if topleftconvention:
            rpc_params['offset_COL'] += 0.5
            rpc_params['offset_LIG'] += 0.5

        return cls(rpc_params)

    @classmethod
    def from_any(cls, primary_file, secondary_file=None, topleftconvention=True):
        """ load from any RPC (auto indetify driver)
        :param primary_file  : rpc filename (dimap, ossim kwl, euclidium coefficients, geotiff)
        :type primary_file  : str
        :param secondary_file  : secondary file (euclidium coefficients)
        :type secondary_file  : str
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        if basename(primary_file).endswith('XML'.upper()):
           dimap_version = identify_dimap(primary_file)
           if dimap_version is not None :
                if float(dimap_version)<2.0 :
                    return cls.from_dimap_v1(primary_file, topleftconvention)
                else:
                    return cls.from_dimap_v2(primary_file, topleftconvention)
        ossim_model = identify_ossim_kwl(primary_file)
        if ossim_model is not None:
            return cls.from_ossim_kwl(primary_file, topleftconvention)
        geotiff_rpc_dict = identify_geotiff_rpc(primary_file)
        if geotiff_rpc_dict is not None:
            return cls.from_geotiff(primary_file, topleftconvention)
        is_eucl_rpc = identify_euclidium_rpc(primary_file)
        if secondary_file is not None:
            is_eucl_rpc = is_eucl_rpc and identify_euclidium_rpc(secondary_file)
        if is_eucl_rpc:
            return cls.from_euclidium(primary_file, secondary_file, topleftconvention)
        ValueError("can''t read rpc file")
        return None

    def calcule_derivees_inv(self,lon,lat,alt):
        """ calcul analytiques des derivees partielles de la loc inverse
            DCdx: derivee de loc_inv_C p/r a X
            DLdy: derivee de loc_inv_L p/r a Y
        """
        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != lon.shape[0]:
            alt = np.full(lon.shape[0], fill_value=alt[0])

        Xnorm = (lon - self.offset_X) / self.scale_X
        Ynorm = (lat - self.offset_Y) / self.scale_Y
        Znorm = (alt - self.offset_ALT) / self.scale_ALT

        if lon.shape[0] > 1:
            DCdx, DCdy, DLdx, DLdy = calcule_derivees_inv_numba(Xnorm, Ynorm, Znorm, np.array(self.Num_COL), np.array(self.Den_COL), np.array(self.Num_LIG), np.array(self.Den_LIG),
                                       self.scale_COL, self.scale_X, self.scale_LIG, self.scale_Y)
        else:
            monomes = self.Monomes[:, 0:1] * np.power(Xnorm, self.Monomes[:, 1:2]) * \
                      np.power(Ynorm, self.Monomes[:, 2:3]) * np.power(Znorm, self.Monomes[:, 3:])

            NumDC = np.dot(self.Num_COL, monomes)
            DenDC = np.dot(self.Den_COL, monomes)
            NumDL = np.dot(self.Num_LIG, monomes)
            DenDL = np.dot(self.Den_LIG, monomes)
            monomes_deriv_x = self.monomes_deriv_1[:, 0:1] * np.power(Xnorm, self.monomes_deriv_1[:, 1:2]) * \
                              np.power(Ynorm, self.monomes_deriv_1[:, 2:3]) * np.power(Znorm, self.monomes_deriv_1[:, 3:])

            NumDCdx = np.dot(self.Num_COL, monomes_deriv_x)
            DenDCdx = np.dot(self.Den_COL, monomes_deriv_x)
            NumDLdx = np.dot(self.Num_LIG, monomes_deriv_x)
            DenDLdx = np.dot(self.Den_LIG, monomes_deriv_x)

            monomes_deriv_y = self.monomes_deriv_2[:, 0:1] * np.power(Xnorm, self.monomes_deriv_2[:, 1:2]) * \
                      np.power(Ynorm, self.monomes_deriv_2[:, 2:3]) * np.power(Znorm, self.monomes_deriv_2[:, 3:])

            NumDCdy = np.dot(self.Num_COL, monomes_deriv_y)
            DenDCdy = np.dot(self.Den_COL, monomes_deriv_y)
            NumDLdy = np.dot(self.Num_LIG, monomes_deriv_y)
            DenDLdy = np.dot(self.Den_LIG, monomes_deriv_y)

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
        print("direct localization not yet implemented for RPC model")
        return None

    def direct_loc_h(self,row,col,alt,fill_nan = False):
        """
        direct localization at constant altitude
        :param row :  line sensor position
        :type row : float or 1D numpy.ndarray dtype=float64
        :param col :  column sensor position
        :type col : float or 1D numpy.ndarray dtype=float64
        :param alt :  altitude
        :param fill_nan : fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan : boolean
        :return ground position (lon,lat,h)
        :rtype numpy.ndarray
        """
        if not isinstance(col, (list, np.ndarray)):
            col = np.array([col])
            row = np.array([row])

        P = np.zeros((col.size, 3))
        filter_nan, P[:,0], P[:,1] = self.filter_coordinates(row, col, fill_nan)
        row = row[filter_nan]
        col = col[filter_nan]

        # Direct localization using direct RPC
        if self.Num_X:
            # ground position
            Xnorm = (col - self.offset_COL)/self.scale_COL
            Ynorm = (row - self.offset_LIG)/self.scale_LIG
            Znorm = (alt - self.offset_ALT)/self.scale_ALT

            if np.sum(abs(Xnorm) > self.lim_extrapol) == Xnorm.shape[0]:
                print("!!!!! l'evaluation au point est extrapolee en colonne ", Xnorm, col)
            if np.sum(abs(Ynorm) > self.lim_extrapol) == Ynorm.shape[0]:
                print("!!!!! l'evaluation au point est extrapolee en ligne ", Ynorm, row)
            if abs(Znorm) > self.lim_extrapol:
                print("!!!!! l'evaluation au point est extrapolee en altitude ", Znorm, alt)

            monomes = self.Monomes[:, 0:1] * np.power(Xnorm, self.Monomes[:, 1:2]) * \
                      np.power(Ynorm, self.Monomes[:, 2:3]) * np.power(Znorm, self.Monomes[:, 3:])

            P[filter_nan, 0] = np.dot(np.array(self.Num_X), monomes)/np.dot(np.array(self.Den_X), monomes)*self.scale_X+self.offset_X
            P[filter_nan, 1] = np.dot(np.array(self.Num_Y), monomes)/np.dot(np.array(self.Den_Y), monomes)*self.scale_Y+self.offset_Y
        # Direct localization using inverse RPC
        else:
            #TODO log info
            print("direct localisation from inverse iterative")
            (P[filter_nan, 0], P[filter_nan, 1], P[filter_nan, 2]) = self.direct_loc_inverse_iterative(row, col, alt, 10, fill_nan)
        P[:, 2] = alt
        return np.squeeze(P)

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
        gri_lon = np.zeros((nbrow,nbcol))
        gri_lat = np.zeros((nbrow,nbcol))
        for c in range(int(nbcol)):
            col = col0 + stepcol*c
            for l in range(int(nbrow)):
                row = row0 + steprow*l
                (gri_lon[l,c],gri_lat[l,c],__) = self.direct_loc_h(row,col,alt)
        return (gri_lon,gri_lat)

    def inverse_loc(self, lon, lat, alt):
        """
        Inverse localization

        :param lon: longitude position
        :type lon : float or 1D numpy.ndarray dtype=float64
        :param lat: latitude position
        :type lat : float or 1D numpy.ndarray dtype=float64
        :param alt: altitude
        :type alt : float
        :return: sensor position (row, col, alt)
        :rtype numpy.ndarray
        """
        if self.Num_COL:
            if not isinstance(lon, (list, np.ndarray)):
                lon = np.array([lon])
                lat = np.array([lat])
                alt = np.array([alt])

            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != lon.shape[0]:
                alt = np.full(lon.shape[0], fill_value=alt[0])

            Xnorm = (lon - self.offset_X) / self.scale_X
            Ynorm = (lat - self.offset_Y) / self.scale_Y
            Znorm = (alt - self.offset_ALT) / self.scale_ALT

            if lon.shape[0] > 1:
                # Inverse localization using numba to reduce calculation time
                Lout, Cout = inverse_loc_numba(Xnorm, Ynorm, Znorm, np.array(self.Num_COL), np.array(self.Den_COL),
                                               np.array(self.Num_LIG), np.array(self.Den_LIG), self.scale_COL,
                                               self.offset_COL, self.scale_LIG, self.offset_LIG)
            else:
                monomes = self.Monomes[:, 0:1] * np.power(Xnorm, self.Monomes[:, 1:2]) * \
                          np.power(Ynorm, self.Monomes[:, 2:3]) * np.power(Znorm, self.Monomes[:, 3:])

                Cout = np.dot(self.Num_COL, monomes) / np.dot(self.Den_COL, monomes) * self.scale_COL + self.offset_COL
                Lout = np.dot(self.Num_LIG, monomes) / np.dot(self.Den_LIG, monomes) * self.scale_LIG + self.offset_LIG
        else:
            print("inverse localisation can't be performed, inverse coefficients have not been defined")
            (Cout, Lout) = (None, None)
        return Lout, Cout, alt

    def filter_coordinates(self, first_coord, second_coord, fill_nan = False, direction = 'direct'):
        """
        Filter nan input values

        :param first_coord :  first coordinate
        :type first_coord : 1D numpy.ndarray dtype=float64
        :param second_coord :  second coordinate
        :type second_coord : 1D numpy.ndarray dtype=float64
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan : boolean
        :param direction :  direct or inverse localisation
        :type direction : str in ('direct', 'inverse')
        :return: filtered coordinates
        :rtype list of numpy.array (index of nan, first filtered, second filtered)
        """
        filter_nan = np.logical_not(np.logical_or(np.isnan(first_coord), np.isnan(second_coord)))

        if fill_nan:
            if direction == 'direct':
                out_x_nan_value = self.offset_X
                out_y_nan_value = self.offset_Y
            else:
                out_x_nan_value = self.offset_COL
                out_y_nan_value = self.offset_LIG
        else:
            out_x_nan_value = np.nan
            out_y_nan_value = np.nan

        x_out = np.full(second_coord.size, out_x_nan_value)
        y_out = np.full(second_coord.size, out_y_nan_value)

        return filter_nan, x_out, y_out

    def direct_loc_inverse_iterative(self, row, col, alt, nb_iter_max=10, fill_nan = False):
        """
        Iterative direct localization using inverse RPC

        :param row :  line sensor position
        :type row : float or 1D numpy.ndarray dtype=float64
        :param col :  column sensor position
        :type col : float or 1D numpy.ndarray dtype=float64
        :param alt :  altitude
        :type alt : float
        :param nb_iter_max: max number of iteration
        :type alt : int
        :param fill_nan: fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan : boolean
        :return: ground position (lon,lat,h)
        :rtype list of numpy.array
        """
        if self.Num_COL:
            if not isinstance(row, (list, np.ndarray)):
                col = np.array([col])
                row = np.array([row])

            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != col.shape[0]:
                alt = np.full(col.shape[0], fill_value=alt[0])

            filter_nan, long_out, lat_out = self.filter_coordinates(row, col, fill_nan)
            row = row[filter_nan]
            col = col[filter_nan]
            alt = alt[filter_nan]

            # if all coord contains Nan then return
            if not np.any(filter_nan):
                return long_out, lat_out, alt

            # inverse localization starting from the center of the scene
            X = np.array([self.offset_X])
            Y = np.array([self.offset_Y])
            (l0, c0, __) = self.inverse_loc(X, Y, alt)

            # desired precision in pixels
            eps = 1e-6

            iteration = 0
            # computing the residue between the sensor positions and those estimated by the inverse localization
            dc = col - c0
            dl = row - l0

            # ground coordinates (latitude and longitude) of each point
            X = np.repeat(X, dc.size)
            Y = np.repeat(Y, dc.size)

            # while the required precision is not achieved
            while (np.max(abs(dc)) > eps or np.max(abs(dl)) > eps) and iteration < nb_iter_max:
                # list of points that require another iteration
                iter_ = np.where((abs(dc) > eps) | (abs(dl) > eps))[0]

                # partial derivatives
                (Cdx, Cdy, Ldx, Ldy) = self.calcule_derivees_inv(X[iter_], Y[iter_], alt[iter_])
                det = Cdx*Ldy-Ldx*Cdy

                dX = (Ldy*dc[iter_] - Cdy*dl[iter_])/det
                dY = (-Ldx*dc[iter_] + Cdx*dl[iter_])/det

                # update ground coordinates
                X[iter_] += dX
                Y[iter_] += dY

                # inverse localization
                (l, c, __) = self.inverse_loc(X[iter_], Y[iter_], alt[iter_])

                # updating the residue between the sensor positions and those estimated by the inverse localization
                dc[iter_] = col[iter_] - c
                dl[iter_] = row[iter_] - l
                iteration += 1

            long_out[filter_nan] = X
            lat_out[filter_nan] = Y

        else:
            print("inverse localisation can't be performed, inverse coefficients have not been defined")
            (long_out, lat_out) = (None, None)

        return long_out, lat_out, alt

    def get_alt_min_max(self):
        """
        returns altitudes min and max layers
        :return alt_min,lat_max
        :rtype list
        """
        return [self.offset_ALT - self.scale_ALT / 2.0, self.offset_ALT + self.scale_ALT / 2.0]

    def los_extrema(self, row, col, alt_min, alt_max, fill_nan = False):
        """
        compute los extrema
        :param row  :  line sensor position
        :type row  : float
        :param col :  column sensor position
        :type col : float
        :param alt_min : los alt min
        :type alt_min  : float
        :param alt_max : los alt max
        :type alt_max : float
        :return los extrema
        :rtype numpy.array (2x3)
        """
        los_edges = np.zeros([2, 3])
        los_edges[0, :] = self.direct_loc_h(row, col, alt_max, fill_nan)
        los_edges[1, :] = self.direct_loc_h(row, col, alt_min, fill_nan)
        return los_edges


@njit('f8(f8, f8, f8, f8[:])', cache=True, fastmath=True)
def polynomial_equation(Xnorm, Ynorm, Znorm, coeff):
    """
    Compute polynomial equation

    :param Xnorm: Normalized longitude position
    :type Xnorm: float 64
    :param Ynorm: Normalized latitude position
    :type Ynorm: float 64
    :param Znorm: Normalized altitude position
    :type Znorm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational
    :rtype: float 64
    """
    rational = coeff[0] + coeff[1] * Xnorm + coeff[2] * Ynorm + coeff[3] * Znorm + coeff[4] * Xnorm * Ynorm + \
        coeff[5] * Xnorm * Znorm + coeff[6] * Ynorm * Znorm + coeff[7] * Xnorm ** 2 + coeff[8] * Ynorm ** 2 + \
        coeff[9] * Znorm ** 2 + coeff[10] * Xnorm * Ynorm * Znorm + coeff[11] * Xnorm ** 3 + \
        coeff[12] * Xnorm * Ynorm ** 2 + coeff[13] * Xnorm * Znorm ** 2 + coeff[14] * Xnorm ** 2 * Ynorm + \
        coeff[15] * Ynorm ** 3 + coeff[16] * Ynorm * Znorm ** 2 + coeff[17] * Xnorm ** 2 * Znorm + \
        coeff[18] * Ynorm ** 2 * Znorm + coeff[19] * Znorm ** 3

    return rational


@njit('Tuple((f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)', parallel=True,
      cache=True, fastmath=True)
def inverse_loc_numba(Xnorm, Ynorm, Znorm, Num_COL, Den_COL, Num_LIG, Den_LIG, scale_COL, offset_COL, scale_LIG,
                  offset_LIG):
    """
    Inverse localization using numba to reduce calculation time on multiple points

    :param Xnorm: Normalized longitude position
    :type Xnorm: 1D np.array dtype np.float 64
    :param Ynorm: Normalized latitude position
    :type Ynorm: 1D np.array dtype np.float 64
    :param Znorm: Normalized altitude position
    :type Znorm: 1D np.array dtype np.float 64
    :param Num_COL: Column numerator coefficients
    :type Num_COL: 1D np.array dtype np.float 64
    :param Den_COL: Column denominator coefficients
    :type Den_COL: 1D np.array dtype np.float 64
    :param Num_LIG: Line numerator coefficients
    :type Num_LIG: 1D np.array dtype np.float 64
    :param Den_LIG: Line denominator coefficients
    :type Den_LIG: 1D np.array dtype np.float 64
    :param scale_COL: Column scale
    :type scale_COL: float 64
    :param offset_COL: Column offset
    :type offset_COL: float 64
    :param scale_LIG: Line scale
    :type scale_LIG: float 64
    :param offset_LIG: Line offset
    :type offset_LIG: float 64
    :return: sensor position (row, col)
    :rtype Tuple(np.ndarray, np.ndarray)
    """
    Cout = np.zeros((Xnorm.shape[0]), dtype=np.float64)
    Lout = np.zeros((Xnorm.shape[0]), dtype=np.float64)

    for i in prange(Xnorm.shape[0]):
        Pu = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Num_COL)
        Qu = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Den_COL)
        Pv = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Num_LIG)
        Qv = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Den_LIG)
        Cout[i] = Pu / Qu * scale_COL + offset_COL
        Lout[i] = Pv / Qv * scale_LIG + offset_LIG

    return Lout, Cout


@njit('f8(f8, f8, f8, f8[:])', cache=True, fastmath=True)
def derivative_polynomial_latitude(Xnorm, Ynorm, Znorm, coeff):
    """
    Compute latitude derivative polynomial equation

    :param Xnorm: Normalized longitude position
    :type Xnorm: float 64
    :param Ynorm: Normalized latitude position
    :type Ynorm: float 64
    :param Znorm: Normalized altitude position
    :type Znorm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    dr = coeff[2] + coeff[4] * Xnorm + coeff[6] * Znorm + 2 * coeff[8] * Ynorm + \
         coeff[10] * Xnorm * Znorm + 2 * coeff[12] * Xnorm * Ynorm + coeff[14] * Xnorm**2 + 3 * coeff[15] * Ynorm**2 + \
         coeff[16] * Znorm**2 + 2 * coeff[18] * Ynorm * Znorm

    return dr


@njit('f8(f8, f8, f8, f8[:])', cache=True, fastmath=True)
def derivative_polynomial_longitude(Xnorm, Ynorm, Znorm, coeff):
    """
    Compute longitude derivative polynomial equation

    :param Xnorm: Normalized longitude position
    :type Xnorm: float 64
    :param Ynorm: Normalized latitude position
    :type Ynorm: float 64
    :param Znorm: Normalized altitude position
    :type Znorm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    dr = coeff[1] + coeff[4] * Ynorm + coeff[5] * Znorm + 2 * coeff[7] * Xnorm + coeff[10] * Ynorm * Znorm + \
         3 * coeff[11] * Xnorm**2 + coeff[12] * Ynorm**2 + coeff[13] * Znorm**2 + 2 * coeff[14] * Ynorm * Xnorm + \
         2 * coeff[17] * Xnorm * Znorm

    return dr


@njit('Tuple((f8[:], f8[:], f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)',
      parallel=True, cache=True, fastmath=True)
def calcule_derivees_inv_numba(Xnorm, Ynorm, Znorm, Num_COL, Den_COL, Num_LIG, Den_LIG, scale_COL, scale_X, scale_LIG, scale_Y):
    """
    Analytically compute the partials derivatives of inverse localization using numba to reduce calculation time on
    multiple points

    :param Xnorm: Normalized longitude position
    :type Xnorm: 1D np.array dtype np.float 64
    :param Ynorm: Normalized latitude position
    :type Ynorm: 1D np.array dtype np.float 64
    :param Znorm: Normalized altitude position
    :type Znorm: 1D np.array dtype np.float 64
    :param Num_COL: Column numerator coefficients
    :type Num_COL: 1D np.array dtype np.float 64
    :param Den_COL: Column denominator coefficients
    :type Den_COL: 1D np.array dtype np.float 64
    :param Num_LIG: Line numerator coefficients
    :type Num_LIG: 1D np.array dtype np.float 64
    :param Den_LIG: Line denominator coefficients
    :type Den_LIG: 1D np.array dtype np.float 64
    :param scale_COL: Column scale
    :type scale_COL: float 64
    :param scale_X: Geodetic longitude scale
    :type scale_X: float 64
    :param scale_LIG: Line scale
    :type scale_LIG: float 64
    :param scale_Y: Geodetic latitude scale
    :type scale_Y: float 64
    :return: partials derivatives of inverse localization
    :rtype: Tuples(DCdx np.array, DCdy np.array, DLdx np.array, DLdy np.array)
    """

    DCdx = np.zeros((Xnorm.shape[0]), dtype=np.float64)
    DCdy = np.zeros((Xnorm.shape[0]), dtype=np.float64)
    DLdx = np.zeros((Xnorm.shape[0]), dtype=np.float64)
    DLdy = np.zeros((Xnorm.shape[0]), dtype=np.float64)

    for i in prange(Xnorm.shape[0]):
        NumDC = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Num_COL)
        DenDC = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Den_COL)
        NumDL = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Num_LIG)
        DenDL = polynomial_equation(Xnorm[i], Ynorm[i], Znorm[i], Den_LIG)

        NumDCdx = derivative_polynomial_longitude(Xnorm[i], Ynorm[i], Znorm[i], Num_COL)
        DenDCdx = derivative_polynomial_longitude(Xnorm[i], Ynorm[i], Znorm[i], Den_COL)
        NumDLdx = derivative_polynomial_longitude(Xnorm[i], Ynorm[i], Znorm[i], Num_LIG)
        DenDLdx = derivative_polynomial_longitude(Xnorm[i], Ynorm[i], Znorm[i], Den_LIG)

        NumDCdy = derivative_polynomial_latitude(Xnorm[i], Ynorm[i], Znorm[i], Num_COL)
        DenDCdy = derivative_polynomial_latitude(Xnorm[i], Ynorm[i], Znorm[i], Den_COL)
        NumDLdy = derivative_polynomial_latitude(Xnorm[i], Ynorm[i], Znorm[i], Num_LIG)
        DenDLdy = derivative_polynomial_latitude(Xnorm[i], Ynorm[i], Znorm[i], Den_LIG)

        DCdx[i] = scale_COL/scale_X * (NumDCdx*DenDC - DenDCdx*NumDC)/DenDC**2
        DCdy[i] = scale_COL/scale_Y * (NumDCdy*DenDC - DenDCdy*NumDC)/DenDC**2
        DLdx[i] = scale_LIG/scale_X * (NumDLdx*DenDL - DenDLdx*NumDL)/DenDL**2
        DLdy[i] = scale_LIG/scale_Y * (NumDLdy*DenDL - DenDLdy*NumDL)/DenDL**2

    return DCdx, DCdy, DLdx, DLdy
