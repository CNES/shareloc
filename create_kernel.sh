export SHARELOCPATH="$( cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P )"
echo "SHARELOCPATH: $SHARELOCPATH"

export CARS_DIR=/work/OT/siaa/3D/Development/guinetj/cars-hal/
cd ${CARS_DIR}
. ./init.sh

cd ${SHARELOCPATH}
pip install -r requirements.txt
pip install -e .

python -m ipykernel install --user --name shareloc --display-name "shareloc_cars"

#copy cars/pandora kernel env in $HOME/.local/share/jupyter/kernels/ 
