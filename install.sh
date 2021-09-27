#!/bin/bash
##
## Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
##
## This file is part of Shareloc
## (see https://github.com/CNES/shareloc).
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

# Set venv name of shareloc venv.
export VENV_NAME="venv"

# Load OTB for gdal and python3
# TODO: move to hpc HAL specific doc or script.
ml otb/7.0-python3.7.2

echo "INFO: clean previous installation..."
make clean
# /tmp noexec configuration.
# TODO: move to hpc HAL specific doc or script.
TMPBUILDDIR=`mktemp -d -p .`
echo "INFO: Shareloc installation..."
TMPDIR=$TMPBUILDDIR VENV=${VENV_NAME} make install
rmdir $TMPBUILDDIR

# Activate virtualenv
source ./${VENV_NAME}/bin/activate
