#!/bin/bash

export SHARELOCPATH="$( cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P )"
echo "SHARELOCPATH: $SHARELOCPATH"

export MODULEPATH=$MODULEPATH:${SHARELOCPATH}/Modulefiles

module load python/3.7.2

# shareloc vars
export PYTHONPATH=$SHARELOCPATH:$PYTHONPATH
export TESTPATH=$SHARELOCPATH/valid_euclidium

# extra application init
export PYENV_HOME=/softs/projets/cars/pyenvs/pyenv-3.7-pandora-v1.b

# virtualenv + autocompletion
source ${PYENV_HOME}/bin/activate
