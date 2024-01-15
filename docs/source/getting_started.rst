.. _getting_started:

===============
Getting started
===============

This section gives a quick easy start on shareloc. 

Please read :ref:`user_manual` for more details, especially :ref:`user_manual_functions`.

Quick Start
===========

* Install Shareloc (section :ref:`install` for more details)

.. code-block:: console

    $ python -m venv shareloc-venv       # Create python virtualenv
    $ source shareloc-venv/bin/activate  # Activate virtualenv
    $ pip install --upgrade pip          # Be sure to have the last PIP version in virtualenv
    $ pip install shareloc               # Install shareloc in virtualenv


* Quick example for direct localization at constant elevation on a grid

1. Get an example of RPC model geom file from shareloc tests

.. code-block:: console
    
    $ wget https://raw.githubusercontent.com/CNES/shareloc/master/tests/data/rpc/phr_ventoux/left_image.geom --no-check-certificate

2. Use shareloc API for direct localization using RPC model
      
.. code-block:: console    

    $ python3
    >>> # Import shareloc modules rpc and localization
    >>> from shareloc.geomodels import GeoModel
    >>> from shareloc.geofunctions.localization import Localization

    >>> # Create RPC object from downloaded geometry file
    >>> rpc_geom_file = "left_image.geom"
    >>> rpc = GeoModel(rpc_geom_file, "RPC") # "RPC" is the geomodel type in ("RPC", "GRID", "RPCoptim") with default value "RPC"

    >>> # Create Localization object from created RPC
    >>> loc = Localization(rpc)

    >>> # Direct localization at first (0, 0) pixel
    >>> loc.direct(0, 0)
    array([[ 5.1608318 , 44.22955181,  0.        ]])

    # --> Result in latitude, longitude, altitude (0 meter over ellipsoid, since altitude is not specified in loc.direct() method)

see :ref:`user_manual_functions` section for more examples.