#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Install CARS, whether via
#      ``python setup.py install``
#    or
#      ``pip install shareloc``

from setuptools import setup
from subprocess import check_output
from codecs import open

# Meta-data.
NAME = 'shareloc'
DESCRIPTION = 'ShareLoc API.'
URL = 'https://gitlab.cnes.fr/OutilsCommuns/shareloc'
AUTHOR = 'CNES'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.1.0'
EMAIL = 'TBD'
LICENSE = 'TBD'
REQUIREMENTS = ['numpy']
DESCRIPTION = '''
ShareLoc API
'''


def readme():
    with open('readme.md', "r", "utf-8") as f:
        return f.read()

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
    python_requires=REQUIRES_PYTHON)






