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

    $ python -m venv venv
    $ source venv/bin/activate
    $ pip install shareloc 


* Quick example for direct localization at constant elevation on a grid

1. Get Grid TIF file from shareloc tests

.. code-block:: console
    
    $ wget https://raw.githubusercontent.com/CNES/shareloc/tests/data/geoide/P1BP--2017030824934340CP/grilles_gld_xH/P1BP--2017030824934340CP.tif    

2. Use shareloc API for direct localization 
      
.. code-block:: console    

    $ python3
    >>> from shareloc.grid import Grid
    >>> from shareloc.localization import Localization
    >>> grid_geom_file = "P1BP--2017030824934340CP.tif"
    >>> grid = Grid(grid_geom_file)
    >>> loc = Localization(grid)
    >>> loc.direct(50.5, 100.5, 100.0)
    array([ 57.21682531,  21.9593944 , 100.        ])
    
    # --> Result in latitude, longitude, altitude.

see :ref:`user_manual_functions` section for more examples.