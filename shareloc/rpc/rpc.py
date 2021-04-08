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
This module contains the RPC class corresponding to the RPC models.
RPC models covered are : DIMAP V1, DIMAP V2, Euclidium, ossim (geom file), geotiff.
"""

from xml.dom import minidom
from os.path import basename
import rasterio as rio
import numpy as np
from numba import njit, prange, config

# Set numba type of threading layer before parallel target compilation
config.THREADING_LAYER = "omp"


def renvoie_linesep(txt_liste_lines):
    """
    Renvoie le separateur de ligne d'un texte sous forme de liste de lignes
    Obtenu par readlines
    """
    if txt_liste_lines[0].endswith("\r\n"):
        line_sep = "\r\n"
    elif txt_liste_lines[0].endswith("\n"):
        line_sep = "\n"
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
    try:
        xmldoc = minidom.parse(xml_file)
        mtd = xmldoc.getElementsByTagName("Metadata_Identification")
        mtd_format = mtd[0].getElementsByTagName("METADATA_FORMAT")[0].firstChild.data
        if mtd_format == "DIMAP_PHR":
            version_tag = "METADATA_PROFILE"
        else:
            version_tag = "METADATA_FORMAT"
        version = mtd[0].getElementsByTagName(version_tag)[0].attributes.items()[0][1]
    except Exception:
        return None
    else:
        return version


def identify_euclidium_rpc(eucl_file):
    """
    parse file to identify if it is an euclidium model (starts with '>>')
    :param eucl_file : euclidium rpc file
    :type eucl_file : str
    :return  True if euclidium rpc has been identified, False otherwise
    :rtype Boolean
    """
    try:
        with open(eucl_file) as euclidium_file:
            content = euclidium_file.readlines()
        is_eucl = True
        for line in content:
            if not line.startswith(">>") and line != "\n":
                is_eucl = False
    except Exception:
        return False
    else:
        return is_eucl


def identify_ossim_kwl(ossim_kwl_file):
    """
    parse geom file to identify if it is an ossim model
    :param ossim_kwl_file : ossim keyword list file
    :type ossim_kwl_file : str
    :return ossim kwl info : ossimmodel or None if not an ossim kwl file
    :rtype str
    """
    try:
        with open(ossim_kwl_file) as ossim_file:
            content = ossim_file.readlines()

        geom_dict = dict()
        for line in content:
            (key, val) = line.split(": ")
            geom_dict[key] = val.rstrip()
        if "type" in geom_dict.keys():
            if geom_dict["type"].strip().startswith("ossim"):
                return geom_dict["type"].strip()
        return None
    except Exception:
        return None


def identify_geotiff_rpc(image_filename):
    """
    read image file to identify if it is a geotiff which contains RPCs
    :param image_filename : image_filename
    :type image_filename : str
    :return rpc info : rpc dict or None  if not a geotiff with rpc
    :rtype str
    """
    try:
        dataset = rio.open(image_filename)
        rpc_dict = dataset.tags(ns="RPC")
        return rpc_dict
    except Exception:
        return None


# gitlab issue #59
# pylint: disable=too-many-branches
def read_eucl_file(eucl_file):
    """
    read euclidium file and parse it
    :param eucl_file : euclidium file
    :type eucl_file : str
    :return parsed file
    :rtype dict
    """
    parsed_file = dict()
    with open(eucl_file, "r") as fid:
        txt = fid.readlines()

    for line in txt:
        if line.startswith(">>\tTYPE_OBJET"):
            if line.split()[-1].endswith("Inverse"):
                parsed_file["type_fic"] = "I"
            if line.split()[-1].endswith("Directe"):
                parsed_file["type_fic"] = "D"

    lsep = renvoie_linesep(txt)

    ind_debut_pxout = txt.index(">>\tCOEFF POLYNOME PXOUT" + lsep)
    ind_debut_qxout = txt.index(">>\tCOEFF POLYNOME QXOUT" + lsep)
    ind_debut_pyout = txt.index(">>\tCOEFF POLYNOME PYOUT" + lsep)
    ind_debut_qyout = txt.index(">>\tCOEFF POLYNOME QYOUT" + lsep)

    coeff_px_str = txt[ind_debut_pxout + 1 : ind_debut_pxout + 21]
    coeff_qx_str = txt[ind_debut_qxout + 1 : ind_debut_qxout + 21]
    coeff_py_str = txt[ind_debut_pyout + 1 : ind_debut_pyout + 21]
    coeff_qy_str = txt[ind_debut_qyout + 1 : ind_debut_qyout + 21]

    poly_coeffs = dict()
    if parsed_file["type_fic"] == "I":
        poly_coeffs["num_col"] = [float(coeff.split()[1]) for coeff in coeff_px_str]
        poly_coeffs["den_col"] = [float(coeff.split()[1]) for coeff in coeff_qx_str]
        poly_coeffs["num_row"] = [float(coeff.split()[1]) for coeff in coeff_py_str]
        poly_coeffs["den_row"] = [float(coeff.split()[1]) for coeff in coeff_qy_str]
    else:
        poly_coeffs["num_x"] = [float(coeff.split()[1]) for coeff in coeff_px_str]
        poly_coeffs["den_x"] = [float(coeff.split()[1]) for coeff in coeff_qx_str]
        poly_coeffs["num_y"] = [float(coeff.split()[1]) for coeff in coeff_py_str]
        poly_coeffs["den_y"] = [float(coeff.split()[1]) for coeff in coeff_qy_str]

    parsed_file["poly_coeffs"] = poly_coeffs
    # list [offset , scale]
    normalisation_coeff = dict()
    for line in txt:
        if line.startswith(">>\tXIN_OFFSET"):
            lsplit = line.split()
            if parsed_file["type_fic"] == "I":
                param = "x"
            else:
                param = "col"
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
        if line.startswith(">>\tYIN_OFFSET"):
            if parsed_file["type_fic"] == "I":
                param = "y"
            else:
                param = "row"
            lsplit = line.split()
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
        if line.startswith(">>\tZIN_OFFSET"):
            lsplit = line.split()
            normalisation_coeff["alt"] = [float(lsplit[4]), float(lsplit[5])]
        if line.startswith(">>\tXOUT_OFFSET"):
            lsplit = line.split()
            if parsed_file["type_fic"] == "D":
                param = "x"
            else:
                param = "col"
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
        if line.startswith(">>\tYOUT_OFFSET"):
            lsplit = line.split()
            if parsed_file["type_fic"] == "D":
                param = "y"
            else:
                param = "row"
            normalisation_coeff[param] = [float(lsplit[4]), float(lsplit[5])]
    parsed_file["normalisation_coeffs"] = normalisation_coeff
    return parsed_file


def check_coeff_consistency(dict1, dict2):
    """
    print an error message inf normalisations coeff are not consistent
    :param dict1 : normalisation coeffs 1
    :type dict1 : dict
    :param dict2 : normalisation coeffs 2
    :type dict2 : dict

    """
    for key, value in dict1.items():
        if dict2[key] != value:
            print(
                "normalisation coeffs are different between"
                " direct en inverse one : {} : {} {}".format(key, value, dict2[key])
            )


class RPC:
    """
    RPC class including direct and inverse localization instance methods
    """

    # gitlab issue #61
    # pylint: disable=too-many-instance-attributes
    def __init__(self, rpc_params):
        for key, value in rpc_params.items():
            setattr(self, key, value)

        self.type = "rpc"
        self.lim_extrapol = 1.0001
        # chaque mononome: c[0]*X**c[1]*Y**c[2]*Z**c[3]
        ordre_monomes_lai = [
            [1, 0, 0, 0],
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [1, 0, 0, 1],
            [1, 1, 1, 0],
            [1, 1, 0, 1],
            [1, 0, 1, 1],
            [1, 2, 0, 0],
            [1, 0, 2, 0],
            [1, 0, 0, 2],
            [1, 1, 1, 1],
            [1, 3, 0, 0],
            [1, 1, 2, 0],
            [1, 1, 0, 2],
            [1, 2, 1, 0],
            [1, 0, 3, 0],
            [1, 0, 1, 2],
            [1, 2, 0, 1],
            [1, 0, 2, 1],
            [1, 0, 0, 3],
        ]

        self.monomes = np.array(ordre_monomes_lai)

        # coefficient des degres monomes avec derivation 1ere variable
        self.monomes_deriv_1 = np.array(
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [1, 0, 1, 0],
                [1, 0, 0, 1],
                [0, 0, 1, 1],
                [2, 1, 0, 0],
                [0, 0, 2, 0],
                [0, 0, 0, 2],
                [1, 0, 1, 1],
                [3, 2, 0, 0],
                [1, 0, 2, 0],
                [1, 0, 0, 2],
                [2, 1, 1, 0],
                [0, 0, 3, 0],
                [0, 0, 1, 2],
                [2, 1, 0, 1],
                [0, 0, 2, 1],
                [0, 0, 0, 3],
            ]
        )

        # coefficient des degres monomes avec derivation 1ere variable
        self.monomes_deriv_2 = np.array(
            [
                [0, 0, 0, 0],
                [0, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 0, 1],
                [1, 1, 0, 0],
                [0, 1, 0, 1],
                [1, 0, 0, 1],
                [0, 2, 0, 0],
                [2, 0, 1, 0],
                [0, 0, 0, 2],
                [1, 1, 0, 1],
                [0, 3, 0, 0],
                [2, 1, 1, 0],
                [0, 1, 0, 2],
                [1, 2, 0, 0],
                [3, 0, 2, 0],
                [1, 0, 0, 2],
                [0, 2, 0, 1],
                [2, 0, 1, 1],
                [0, 0, 0, 3],
            ]
        )

        self.inverse_coefficient = False
        self.direct_coefficient = False

        # pylint: disable=access-member-before-definition
        if self.num_col:
            self.inverse_coefficient = True
            self.num_col = np.array(self.num_col)
            self.den_col = np.array(self.den_col)
            self.num_row = np.array(self.num_row)
            self.den_row = np.array(self.den_row)

        # pylint: disable=access-member-before-definition
        if self.num_x:
            self.direct_coefficient = True
            self.num_x = np.array(self.num_x)
            self.den_x = np.array(self.den_x)
            self.num_y = np.array(self.num_y)
            self.den_y = np.array(self.den_y)

        self.alt_minmax = [self.offset_alt - self.scale_alt, self.offset_alt + self.scale_alt]

    @classmethod
    def from_dimap(cls, dimap_filepath, topleftconvention=True):
        """load from dimap
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
            if float(dimap_version) >= 2.0:
                return cls.from_dimap_v2(dimap_filepath, topleftconvention)
        else:
            ValueError("can''t read dimap file")
            return None

    @classmethod
    def from_dimap_v2(cls, dimap_filepath, topleftconvention=True):
        """load from dimap  v2
        :param dimap_filepath  : dimap xml file
        :type dimap_filepath  : str
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """

        rpc_params = dict()

        if not basename(dimap_filepath).endswith("XML".upper()):
            raise ValueError("dimap must ends with .xml")

        xmldoc = minidom.parse(dimap_filepath)

        mtd = xmldoc.getElementsByTagName("Metadata_Identification")
        version = mtd[0].getElementsByTagName("METADATA_FORMAT")[0].attributes.items()[0][1]
        rpc_params["driver_type"] = "dimap_v" + version
        global_rfm = xmldoc.getElementsByTagName("Global_RFM")[0]
        direct_coeffs = global_rfm.getElementsByTagName("Direct_Model")[0]
        rpc_params["num_x"] = [
            float(direct_coeffs.getElementsByTagName("SAMP_NUM_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_x"] = [
            float(direct_coeffs.getElementsByTagName("SAMP_DEN_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["num_y"] = [
            float(direct_coeffs.getElementsByTagName("LINE_NUM_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_y"] = [
            float(direct_coeffs.getElementsByTagName("LINE_DEN_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        inverse_coeffs = global_rfm.getElementsByTagName("Inverse_Model")[0]
        rpc_params["num_col"] = [
            float(inverse_coeffs.getElementsByTagName("SAMP_NUM_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_col"] = [
            float(inverse_coeffs.getElementsByTagName("SAMP_DEN_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["num_row"] = [
            float(inverse_coeffs.getElementsByTagName("LINE_NUM_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        rpc_params["den_row"] = [
            float(inverse_coeffs.getElementsByTagName("LINE_DEN_COEFF_{}".format(index))[0].firstChild.data)
            for index in range(1, 21)
        ]
        normalisation_coeffs = global_rfm.getElementsByTagName("RFM_Validity")[0]
        rpc_params["offset_col"] = float(normalisation_coeffs.getElementsByTagName("SAMP_OFF")[0].firstChild.data)
        rpc_params["scale_col"] = float(normalisation_coeffs.getElementsByTagName("SAMP_SCALE")[0].firstChild.data)
        rpc_params["offset_row"] = float(normalisation_coeffs.getElementsByTagName("LINE_OFF")[0].firstChild.data)
        rpc_params["scale_row"] = float(normalisation_coeffs.getElementsByTagName("LINE_SCALE")[0].firstChild.data)
        rpc_params["offset_alt"] = float(normalisation_coeffs.getElementsByTagName("HEIGHT_OFF")[0].firstChild.data)
        rpc_params["scale_alt"] = float(normalisation_coeffs.getElementsByTagName("HEIGHT_SCALE")[0].firstChild.data)
        rpc_params["offset_x"] = float(normalisation_coeffs.getElementsByTagName("LONG_OFF")[0].firstChild.data)
        rpc_params["scale_x"] = float(normalisation_coeffs.getElementsByTagName("LONG_SCALE")[0].firstChild.data)
        rpc_params["offset_y"] = float(normalisation_coeffs.getElementsByTagName("LAT_OFF")[0].firstChild.data)
        rpc_params["scale_y"] = float(normalisation_coeffs.getElementsByTagName("LAT_SCALE")[0].firstChild.data)
        rpc_params["offset_col"] -= 1.0
        rpc_params["offset_row"] -= 1.0
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params["offset_col"] += 0.5
            rpc_params["offset_row"] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_dimap_v1(cls, dimap_filepath, topleftconvention=True):
        """load from dimap  v1
        :param dimap_filepath  : dimap xml file
        :type dimap_filepath  : str
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        rpc_params = dict()

        if not basename(dimap_filepath).endswith("XML".upper()):
            raise ValueError("dimap must ends with .xml")

        xmldoc = minidom.parse(dimap_filepath)

        mtd = xmldoc.getElementsByTagName("Metadata_Identification")
        version = mtd[0].getElementsByTagName("METADATA_PROFILE")[0].attributes.items()[0][1]
        rpc_params["driver_type"] = "dimap_v" + version

        global_rfm = xmldoc.getElementsByTagName("Global_RFM")
        rfm_validity = xmldoc.getElementsByTagName("RFM_Validity")
        coeff_lon = [float(el) for el in global_rfm[0].getElementsByTagName("F_LON")[0].firstChild.data.split()]
        coeff_lat = [float(el) for el in global_rfm[0].getElementsByTagName("F_LAT")[0].firstChild.data.split()]
        coeff_col = [float(el) for el in global_rfm[0].getElementsByTagName("F_COL")[0].firstChild.data.split()]
        coeff_lig = [float(el) for el in global_rfm[0].getElementsByTagName("F_ROW")[0].firstChild.data.split()]

        scale_lon = float(rfm_validity[0].getElementsByTagName("Lon")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_lon = float(rfm_validity[0].getElementsByTagName("Lon")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_lat = float(rfm_validity[0].getElementsByTagName("Lat")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_lat = float(rfm_validity[0].getElementsByTagName("Lat")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_alt = float(rfm_validity[0].getElementsByTagName("Alt")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_alt = float(rfm_validity[0].getElementsByTagName("Alt")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_col = float(rfm_validity[0].getElementsByTagName("Col")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_col = float(rfm_validity[0].getElementsByTagName("Col")[0].getElementsByTagName("B")[0].firstChild.data)
        scale_row = float(rfm_validity[0].getElementsByTagName("Row")[0].getElementsByTagName("A")[0].firstChild.data)
        offset_row = float(rfm_validity[0].getElementsByTagName("Row")[0].getElementsByTagName("B")[0].firstChild.data)

        rpc_params["offset_col"] = offset_col
        rpc_params["scale_col"] = scale_col
        rpc_params["offset_row"] = offset_row
        rpc_params["scale_row"] = scale_row
        rpc_params["offset_alt"] = offset_alt
        rpc_params["scale_alt"] = scale_alt
        rpc_params["offset_x"] = offset_lon
        rpc_params["scale_x"] = scale_lon
        rpc_params["offset_y"] = offset_lat
        rpc_params["scale_y"] = scale_lat
        rpc_params["num_x"] = coeff_lon[0:20]
        rpc_params["den_x"] = coeff_lon[20::]
        rpc_params["num_y"] = coeff_lat[0:20]
        rpc_params["den_y"] = coeff_lat[20::]
        rpc_params["num_col"] = coeff_col[0:20]
        rpc_params["den_col"] = coeff_col[20::]
        rpc_params["num_row"] = coeff_lig[0:20]
        rpc_params["den_row"] = coeff_lig[20::]
        rpc_params["offset_col"] -= 1.0
        rpc_params["offset_row"] -= 1.0
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params["offset_col"] += 0.5
            rpc_params["offset_row"] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_geotiff(cls, image_filename, topleftconvention=True):
        """Load from a  geotiff image file
        :param image_filename  : image filename
        :type image_filename  : str
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        dataset = rio.open(image_filename)
        rpc_dict = dataset.tags(ns="RPC")
        if not rpc_dict:
            print("{} does not contains RPCs ".format(image_filename))
            raise ValueError
        rpc_params = dict()
        rpc_params["den_row"] = parse_coeff_line(rpc_dict["LINE_DEN_COEFF"])
        rpc_params["num_row"] = parse_coeff_line(rpc_dict["LINE_NUM_COEFF"])
        rpc_params["num_col"] = parse_coeff_line(rpc_dict["SAMP_NUM_COEFF"])
        rpc_params["den_col"] = parse_coeff_line(rpc_dict["SAMP_DEN_COEFF"])
        rpc_params["offset_col"] = float(rpc_dict["SAMP_OFF"])
        rpc_params["scale_col"] = float(rpc_dict["SAMP_SCALE"])
        rpc_params["offset_row"] = float(rpc_dict["LINE_OFF"])
        rpc_params["scale_row"] = float(rpc_dict["LINE_SCALE"])
        rpc_params["offset_alt"] = float(rpc_dict["HEIGHT_OFF"])
        rpc_params["scale_alt"] = float(rpc_dict["HEIGHT_SCALE"])
        rpc_params["offset_x"] = float(rpc_dict["LONG_OFF"])
        rpc_params["scale_x"] = float(rpc_dict["LONG_SCALE"])
        rpc_params["offset_y"] = float(rpc_dict["LAT_OFF"])
        rpc_params["scale_y"] = float(rpc_dict["LAT_SCALE"])
        # inverse coeff are not defined
        rpc_params["num_x"] = None
        rpc_params["den_x"] = None
        rpc_params["num_y"] = None
        rpc_params["den_y"] = None
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params["offset_col"] += 0.5
            rpc_params["offset_row"] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_ossim_kwl(cls, ossim_kwl_filename, topleftconvention=True):
        """Load from a geom file
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        rpc_params = dict()
        # OSSIM keyword list
        rpc_params["driver_type"] = "ossim_kwl"
        with open(ossim_kwl_filename) as ossim_file:
            content = ossim_file.readlines()

        geom_dict = dict()
        for line in content:
            (key, val) = line.split(": ")
            geom_dict[key] = val.rstrip()

        rpc_params["den_row"] = [np.nan] * 20
        rpc_params["num_row"] = [np.nan] * 20
        rpc_params["den_col"] = [np.nan] * 20
        rpc_params["num_col"] = [np.nan] * 20
        for index in range(0, 20):
            axis = "line"
            num_den = "den"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params["den_row"][index] = float(geom_dict[key])
            num_den = "num"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params["num_row"][index] = float(geom_dict[key])
            axis = "samp"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params["num_col"][index] = float(geom_dict[key])
            num_den = "den"
            key = "{0}_{1}_coeff_{2:02d}".format(axis, num_den, index)
            rpc_params["den_col"][index] = float(geom_dict[key])
        rpc_params["offset_col"] = float(geom_dict["samp_off"])
        rpc_params["scale_col"] = float(geom_dict["samp_scale"])
        rpc_params["offset_row"] = float(geom_dict["line_off"])
        rpc_params["scale_row"] = float(geom_dict["line_scale"])
        rpc_params["offset_alt"] = float(geom_dict["height_off"])
        rpc_params["scale_alt"] = float(geom_dict["height_scale"])
        rpc_params["offset_x"] = float(geom_dict["long_off"])
        rpc_params["scale_x"] = float(geom_dict["long_scale"])
        rpc_params["offset_y"] = float(geom_dict["lat_off"])
        rpc_params["scale_y"] = float(geom_dict["lat_scale"])
        # inverse coeff are not defined
        rpc_params["num_x"] = None
        rpc_params["den_x"] = None
        rpc_params["num_y"] = None
        rpc_params["den_y"] = None
        # If top left convention, 0.5 pixel shift added on col/row offsets
        if topleftconvention:
            rpc_params["offset_col"] += 0.5
            rpc_params["offset_row"] += 0.5
        return cls(rpc_params)

    @classmethod
    def from_euclidium(cls, primary_euclidium_coeff, secondary_euclidium_coeff=None, topleftconvention=True):
        """load from euclidium
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
        rpc_params["driver_type"] = "euclidium"

        # lecture fichier euclidium
        primary_coeffs = read_eucl_file(primary_euclidium_coeff)

        # info log
        print("primary euclidium file is of {} type".format(primary_coeffs["type_fic"]))

        rpc_params["num_x"] = None
        rpc_params["den_x"] = None
        rpc_params["num_y"] = None
        rpc_params["den_y"] = None
        rpc_params["num_col"] = None
        rpc_params["den_col"] = None
        rpc_params["num_row"] = None
        rpc_params["den_row"] = None

        for key, value in primary_coeffs["poly_coeffs"].items():
            rpc_params[key] = value

        rpc_params["normalisation_coeffs"] = primary_coeffs["normalisation_coeffs"]
        for key, value in primary_coeffs["normalisation_coeffs"].items():
            rpc_params["offset_" + key] = value[0]
            rpc_params["scale_" + key] = value[1]

        if secondary_euclidium_coeff is not None:
            secondary_coeffs = read_eucl_file(secondary_euclidium_coeff)
            print("secondary euclidium file is of {} type".format(secondary_coeffs["type_fic"]))
            check_coeff_consistency(primary_coeffs["normalisation_coeffs"], secondary_coeffs["normalisation_coeffs"])

            for key, value in secondary_coeffs["poly_coeffs"].items():
                rpc_params[key] = value

        rpc_params["offset_col"] -= 0.5
        rpc_params["offset_row"] -= 0.5
        # If top left convention, 0.5 pixel shift added on col/row offsets

        if topleftconvention:
            rpc_params["offset_col"] += 0.5
            rpc_params["offset_row"] += 0.5

        return cls(rpc_params)

    @classmethod
    def from_any(cls, primary_file, secondary_file=None, topleftconvention=True):
        """load from any RPC (auto indetify driver)
        :param primary_file  : rpc filename (dimap, ossim kwl, euclidium coefficients, geotiff)
        :type primary_file  : str
        :param secondary_file  : secondary file (euclidium coefficients)
        :type secondary_file  : str
        :param topleftconvention  : [0,0] position
        :type topleftconvention  : boolean
        If False : [0,0] is at the center of the Top Left pixel
        If True : [0,0] is at the top left of the Top Left pixel (OSSIM)
        """
        if basename(primary_file).endswith("XML".upper()):
            dimap_version = identify_dimap(primary_file)
            if dimap_version is not None:
                if float(dimap_version) < 2.0:
                    return cls.from_dimap_v1(primary_file, topleftconvention)
                if float(dimap_version) >= 2.0:
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
        ValueError("can not read rpc file")
        return None

    def calcule_derivees_inv(self, lon, lat, alt):
        """calcul analytiques des derivees partielles de la loc inverse
        DCdx: derivee de loc_inv_C p/r a X
        DLdy: derivee de loc_inv_L p/r a Y
        """
        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != lon.shape[0]:
            alt = np.full(lon.shape[0], fill_value=alt[0])

        lon_norm = (lon - self.offset_x) / self.scale_x
        lat_norm = (lat - self.offset_y) / self.scale_y
        alt_norm = (alt - self.offset_alt) / self.scale_alt

        dcol_dlon, dcol_dlat, drow_dlon, drow_dlat = calcule_derivees_inv_numba(
            lon_norm,
            lat_norm,
            alt_norm,
            self.num_col,
            self.den_col,
            self.num_row,
            self.den_row,
            self.scale_col,
            self.scale_x,
            self.scale_row,
            self.scale_y,
        )

        return (dcol_dlon, dcol_dlat, drow_dlon, drow_dlat)

    def direct_loc_h(self, row, col, alt, fill_nan=False):
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

        if not isinstance(alt, (list, np.ndarray)):
            alt = np.array([alt])

        if alt.shape[0] != col.shape[0]:
            alt = np.full(col.shape[0], fill_value=alt[0])

        points = np.zeros((col.size, 3))
        filter_nan, points[:, 0], points[:, 1] = self.filter_coordinates(row, col, fill_nan)
        row = row[filter_nan]
        col = col[filter_nan]

        # Direct localization using direct RPC
        if self.direct_coefficient:
            # ground position
            col_norm = (col - self.offset_col) / self.scale_col
            row_norm = (row - self.offset_row) / self.scale_row
            alt_norm = (alt - self.offset_alt) / self.scale_alt

            points[filter_nan, 1], points[filter_nan, 0] = compute_rational_function_polynomial(
                col_norm,
                row_norm,
                alt_norm,
                self.num_x,
                self.den_x,
                self.num_y,
                self.den_y,
                self.scale_x,
                self.offset_x,
                self.scale_y,
                self.offset_y,
            )

        # Direct localization using inverse RPC
        else:
            # print("direct localisation from inverse iterative")
            (points[filter_nan, 0], points[filter_nan, 1], points[filter_nan, 2]) = self.direct_loc_inverse_iterative(
                row, col, alt, 10, fill_nan
            )
        points[:, 2] = alt
        return np.squeeze(points)

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
        gri_lon = np.zeros((nbrow, nbcol))
        gri_lat = np.zeros((nbrow, nbcol))
        for column in range(int(nbcol)):
            col = col0 + stepcol * column
            for line in range(int(nbrow)):
                row = row0 + steprow * line
                (gri_lon[line, column], gri_lat[line, column], __) = self.direct_loc_h(row, col, alt)
        return (gri_lon, gri_lat)

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
        if isinstance(col, (list, np.ndarray)):
            points_nb = len(col)
        else:
            points_nb = 1
            row = np.array([row])
            col = np.array([col])
        direct_dtm = np.zeros((points_nb, 3))

        for i in range(points_nb):
            row_i = row[i]
            col_i = col[i]
            # print("min {} max {}".format(dtm.Zmin,dtm.Zmax))
            (min_dtm, max_dtm) = (dtm.alt_min - 1.0, dtm.alt_max + 1.0)
            # if min_dtm < self.offset_alt - self.scale_alt:
            #    print("minimum dtm value is outside RPC validity domain")
            # if max_dtm > self.offset_alt + self.scale_alt:
            #    print("maximum dtm value is outside RPC validity domain")
            los = self.los_extrema(row_i, col_i, min_dtm, max_dtm)
            (__, __, position_cube, alti) = dtm.intersect_dtm_cube(los)
            (__, __, position) = dtm.intersection(los, position_cube, alti)
            direct_dtm[i, :] = position
        return direct_dtm

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
        if self.inverse_coefficient:
            if not isinstance(lon, (list, np.ndarray)):
                lon = np.array([lon])
                lat = np.array([lat])
            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if not isinstance(alt, (list, np.ndarray)):
                alt = np.array([alt])

            if alt.shape[0] != lon.shape[0]:
                alt = np.full(lon.shape[0], fill_value=alt[0])

            lon_norm = (lon - self.offset_x) / self.scale_x
            lat_norm = (lat - self.offset_y) / self.scale_y
            alt_norm = (alt - self.offset_alt) / self.scale_alt

            row_out, col_out = compute_rational_function_polynomial(
                lon_norm,
                lat_norm,
                alt_norm,
                self.num_col,
                self.den_col,
                self.num_row,
                self.den_row,
                self.scale_col,
                self.offset_col,
                self.scale_row,
                self.offset_row,
            )
        else:
            print("inverse localisation can't be performed, inverse coefficients have not been defined")
            (col_out, row_out) = (None, None)
        return row_out, col_out, alt

    def filter_coordinates(self, first_coord, second_coord, fill_nan=False, direction="direct"):
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
            if direction == "direct":
                out_x_nan_value = self.offset_x
                out_y_nan_value = self.offset_y
            else:
                out_x_nan_value = self.offset_col
                out_y_nan_value = self.offset_row
        else:
            out_x_nan_value = np.nan
            out_y_nan_value = np.nan

        x_out = np.full(len(second_coord), out_x_nan_value)
        y_out = np.full(len(second_coord), out_y_nan_value)

        return filter_nan, x_out, y_out

    def direct_loc_inverse_iterative(self, row, col, alt, nb_iter_max=10, fill_nan=False):
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
        if self.inverse_coefficient:
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
            lon = np.array([self.offset_x])
            lat = np.array([self.offset_y])
            (row_start, col_start, __) = self.inverse_loc(lon, lat, alt)

            # desired precision in pixels
            eps = 1e-6

            iteration = 0
            # computing the residue between the sensor positions and those estimated by the inverse localization
            delta_col = col - col_start
            delta_row = row - row_start

            # ground coordinates (latitude and longitude) of each point
            lon = np.repeat(lon, delta_col.size)
            lat = np.repeat(lat, delta_col.size)
            # while the required precision is not achieved
            while (np.max(abs(delta_col)) > eps or np.max(abs(delta_row)) > eps) and iteration < nb_iter_max:
                # list of points that require another iteration
                iter_ = np.where((abs(delta_col) > eps) | (abs(delta_row) > eps))[0]

                # partial derivatives
                (dcol_dlon, dcol_dlat, drow_dlon, drow_dlat) = self.calcule_derivees_inv(
                    lon[iter_], lat[iter_], alt[iter_]
                )
                det = dcol_dlon * drow_dlat - drow_dlon * dcol_dlat

                delta_lon = (drow_dlat * delta_col[iter_] - dcol_dlat * delta_row[iter_]) / det
                delta_lat = (-drow_dlon * delta_col[iter_] + dcol_dlon * delta_row[iter_]) / det

                # update ground coordinates
                lon[iter_] += delta_lon
                lat[iter_] += delta_lat

                # inverse localization
                (row_estim, col_estim, __) = self.inverse_loc(lon[iter_], lat[iter_], alt[iter_])

                # updating the residue between the sensor positions and those estimated by the inverse localization
                delta_col[iter_] = col[iter_] - col_estim
                delta_row[iter_] = row[iter_] - row_estim
                iteration += 1

            long_out[filter_nan] = lon
            lat_out[filter_nan] = lat

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
        return [self.offset_alt - self.scale_alt / 2.0, self.offset_alt + self.scale_alt / 2.0]

    def los_extrema(self, row, col, alt_min=None, alt_max=None, fill_nan=False):
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
        extrapolate = False
        if alt_min is None or alt_max is None:
            [los_alt_min, los_alt_max] = self.get_alt_min_max()
        elif alt_min >= self.alt_minmax[0] and alt_max <= self.alt_minmax[1]:
            los_alt_min = alt_min
            los_alt_max = alt_max
        else:
            extrapolate = True
            [los_alt_min, los_alt_max] = self.get_alt_min_max()

        los_edges = np.zeros([2, 3])
        los_edges = self.direct_loc_h(
            np.array([row, row]), np.array([col, col]), np.array([los_alt_max, los_alt_min]), fill_nan
        )
        if extrapolate:
            diff = los_edges[0, :] - los_edges[1, :]
            delta_alt = diff[2]
            los_edges[0, :] = los_edges[1, :] + diff * (alt_max - los_edges[1, 2]) / delta_alt
            los_edges[1, :] = los_edges[1, :] + diff * (alt_min - los_edges[1, 2]) / delta_alt

        return los_edges


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def polynomial_equation(xnorm, ynorm, znorm, coeff):
    """
    Compute polynomial equation

    :param xnorm: Normalized longitude (for inverse) or column (for direct) position
    :type xnorm: float 64
    :param ynorm: Normalized latitude (for inverse) or line (for direct) position
    :type ynorm: float 64
    :param znorm: Normalized altitude position
    :type znorm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational
    :rtype: float 64
    """
    rational = (
        coeff[0]
        + coeff[1] * xnorm
        + coeff[2] * ynorm
        + coeff[3] * znorm
        + coeff[4] * xnorm * ynorm
        + coeff[5] * xnorm * znorm
        + coeff[6] * ynorm * znorm
        + coeff[7] * xnorm ** 2
        + coeff[8] * ynorm ** 2
        + coeff[9] * znorm ** 2
        + coeff[10] * xnorm * ynorm * znorm
        + coeff[11] * xnorm ** 3
        + coeff[12] * xnorm * ynorm ** 2
        + coeff[13] * xnorm * znorm ** 2
        + coeff[14] * xnorm ** 2 * ynorm
        + coeff[15] * ynorm ** 3
        + coeff[16] * ynorm * znorm ** 2
        + coeff[17] * xnorm ** 2 * znorm
        + coeff[18] * ynorm ** 2 * znorm
        + coeff[19] * znorm ** 3
    )

    return rational


# pylint: disable=too-many-arguments
@njit(
    "Tuple((f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)",
    parallel=True,
    cache=True,
    fastmath=True,
)
def compute_rational_function_polynomial(
    lon_col_norm,
    lat_row_norm,
    alt_norm,
    num_col,
    den_col,
    num_lin,
    den_lin,
    scale_col,
    offset_col,
    scale_lin,
    offset_lin,
):
    """
    Compute rational function polynomial using numba to reduce calculation time on multiple points.
    useful to compute direct and inverse localization using direct or inverse RPC.

    :param lon_col_norm: Normalized longitude (for inverse) or column (for direct) position
    :type lon_col_norm: 1D np.array dtype np.float 64
    :param lat_row_norm: Normalized latitude (for inverse) or line (for direct) position
    :type lat_row_norm: 1D np.array dtype np.float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: 1D np.array dtype np.float 64
    :param num_col: Column numerator coefficients
    :type num_col: 1D np.array dtype np.float 64
    :param den_col: Column denominator coefficients
    :type den_col: 1D np.array dtype np.float 64
    :param num_lin: Line numerator coefficients
    :type num_lin: 1D np.array dtype np.float 64
    :param den_lin: Line denominator coefficients
    :type den_lin: 1D np.array dtype np.float 64
    :param scale_col: Column scale
    :type scale_col: float 64
    :param offset_col: Column offset
    :type offset_col: float 64
    :param scale_lin: Line scale
    :type scale_lin: float 64
    :param offset_lin: Line offset
    :type offset_lin: float 64
    :return: for inverse localization : sensor position (row, col). for direct localization : ground position (lon, lat)
    :rtype Tuple(np.ndarray, np.ndarray)
    """
    assert lon_col_norm.shape == alt_norm.shape

    col_lat_out = np.zeros((lon_col_norm.shape[0]), dtype=np.float64)
    row_lon_out = np.zeros((lon_col_norm.shape[0]), dtype=np.float64)

    # pylint: disable=not-an-iterable
    for i in prange(lon_col_norm.shape[0]):
        poly_num_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_col)
        poly_den_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_col)
        poly_num_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_lin)
        poly_den_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_lin)
        col_lat_out[i] = poly_num_col / poly_den_col * scale_col + offset_col
        row_lon_out[i] = poly_num_lin / poly_den_lin * scale_lin + offset_lin

    return row_lon_out, col_lat_out


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, coeff):
    """
    Compute latitude derivative polynomial equation

    :param lon_norm: Normalized longitude position
    :type lon_norm: float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    derivate = (
        coeff[2]
        + coeff[4] * lon_norm
        + coeff[6] * alt_norm
        + 2 * coeff[8] * lat_norm
        + coeff[10] * lon_norm * alt_norm
        + 2 * coeff[12] * lon_norm * lat_norm
        + coeff[14] * lon_norm ** 2
        + 3 * coeff[15] * lat_norm ** 2
        + coeff[16] * alt_norm ** 2
        + 2 * coeff[18] * lat_norm * alt_norm
    )

    return derivate


@njit("f8(f8, f8, f8, f8[:])", cache=True, fastmath=True)
def derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, coeff):
    """
    Compute longitude derivative polynomial equation

    :param lon_norm: Normalized longitude position
    :type lon_norm: float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: float 64
    :param coeff: coefficients
    :type coeff: 1D np.array dtype np.float 64
    :return: rational derivative
    :rtype: float 64
    """
    derivate = (
        coeff[1]
        + coeff[4] * lat_norm
        + coeff[5] * alt_norm
        + 2 * coeff[7] * lon_norm
        + coeff[10] * lat_norm * alt_norm
        + 3 * coeff[11] * lon_norm ** 2
        + coeff[12] * lat_norm ** 2
        + coeff[13] * alt_norm ** 2
        + 2 * coeff[14] * lat_norm * lon_norm
        + 2 * coeff[17] * lon_norm * alt_norm
    )

    return derivate


# pylint: disable=too-many-arguments
@njit(
    "Tuple((f8[:], f8[:], f8[:], f8[:]))(f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8, f8, f8, f8)",
    parallel=True,
    cache=True,
    fastmath=True,
)
def calcule_derivees_inv_numba(
    lon_norm, lat_norm, alt_norm, num_col, den_col, num_lin, den_lin, scale_col, scale_lon, scale_lin, scale_lat
):
    """
    Analytically compute the partials derivatives of inverse localization using numba to reduce calculation time on
    multiple points

    :param lon_norm: Normalized longitude position
    :type lon_norm: 1D np.array dtype np.float 64
    :param lat_norm: Normalized latitude position
    :type lat_norm: 1D np.array dtype np.float 64
    :param alt_norm: Normalized altitude position
    :type alt_norm: 1D np.array dtype np.float 64
    :param num_col: Column numerator coefficients
    :type num_col: 1D np.array dtype np.float 64
    :param den_col: Column denominator coefficients
    :type den_col: 1D np.array dtype np.float 64
    :param num_lin: Line numerator coefficients
    :type num_lin: 1D np.array dtype np.float 64
    :param den_lin: Line denominator coefficients
    :type den_lin: 1D np.array dtype np.float 64
    :param scale_col: Column scale
    :type scale_col: float 64
    :param scale_lon: Geodetic longitude scale
    :type scale_lon: float 64
    :param scale_lin: Line scale
    :type scale_lin: float 64
    :param scale_lat: Geodetic latitude scale
    :type scale_lat: float 64
    :return: partials derivatives of inverse localization
    :rtype: Tuples(dcol_dlon np.array, dcol_dlat np.array, drow_dlon np.array, drow_dlat np.array)
    """

    dcol_dlon = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    dcol_dlat = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    drow_dlon = np.zeros((lon_norm.shape[0]), dtype=np.float64)
    drow_dlat = np.zeros((lon_norm.shape[0]), dtype=np.float64)

    # pylint: disable=not-an-iterable
    for i in prange(lon_norm.shape[0]):
        num_dcol = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow = polynomial_equation(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        num_dcol_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow_dlon = derivative_polynomial_longitude(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        num_dcol_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], num_col)
        den_dcol_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], den_col)
        num_drow_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], num_lin)
        den_drow_dlat = derivative_polynomial_latitude(lon_norm[i], lat_norm[i], alt_norm[i], den_lin)

        dcol_dlon[i] = scale_col / scale_lon * (num_dcol_dlon * den_dcol - den_dcol_dlon * num_dcol) / den_dcol ** 2
        dcol_dlat[i] = scale_col / scale_lat * (num_dcol_dlat * den_dcol - den_dcol_dlat * num_dcol) / den_dcol ** 2
        drow_dlon[i] = scale_lin / scale_lon * (num_drow_dlon * den_drow - den_drow_dlon * num_drow) / den_drow ** 2
        drow_dlat[i] = scale_lin / scale_lat * (num_drow_dlat * den_drow - den_drow_dlat * num_drow) / den_drow ** 2

    return dcol_dlon, dcol_dlat, drow_dlon, drow_dlat
