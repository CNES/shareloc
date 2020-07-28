#!/bin/bash

export SHARELOCPATH="$( cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P )"
echo "SHARELOCPATH: $SHARELOCPATH"

export TESTPATH=$SHARELOCPATH/valid

ml python/3.7.2
virtualenv venv_shareloc
source ./venv_shareloc/bin/activate
pip install -e .
pip install pytest
