#!/bin/bash

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
