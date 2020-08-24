#!/bin/bash

export SHARELOCPATH="$( cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P )"
echo "SHARELOCPATH: $SHARELOCPATH"

export TESTPATH=$SHARELOCPATH/valid

PYTHON_VERSION=3.7.2
SHARELOC_VERSION=0.1.0
VENV_NAME=shareloc-${SHARELOC_VERSION}-pyenv-${PYTHON_VERSION}
ml python/$PYTHON_VERSION
virtualenv ${VENV_NAME}
source ./${VENV_NAME}/bin/activate
pip install -e .
pip install pytest
