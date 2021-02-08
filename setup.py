#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Install CARS, whether via
#      ``python setup.py install``
#    or
#      ``pip install shareloc``
"""
This module contains the required libraries and softwares allowing to execute the software,
and setup elements to configure and identify the software.
"""

from codecs import open as copen
from setuptools import setup, find_packages

# Meta-data.
NAME = 'shareloc'
DESCRIPTION = 'ShareLoc API.'
URL = 'https://gitlab.cnes.fr/OutilsCommuns/shareloc'
AUTHOR = 'CNES'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.1.0'
EMAIL = 'TBD'
LICENSE = 'TBD'
REQUIREMENTS = ['numpy','gdal','rasterio','xarray','netCDF4', 'numba']
DESCRIPTION = '''
ShareLoc API
'''


def readme():
    with copen('readme.md', 'r', 'utf-8') as fstream:
        return fstream.read()

# Setup
setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    url=URL,
    author=AUTHOR,
    author_email=EMAIL,
    license=LICENSE,
    long_description=DESCRIPTION,
    install_requires=REQUIREMENTS,
    python_requires=REQUIRES_PYTHON,
    packages=find_packages())






