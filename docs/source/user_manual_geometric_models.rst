.. _user_manual_geometric_models:


================
Geometric models
================

Shareloc handles two type of geometric models, RPC and Direct location grids.

RPC
===

:term:`RPC` is an analytics function for ground (lon,lat,h) to image (r,c) mapping, it can be summarized as :math:`(r,c) = (f(P,L,H),g(P,L,H))`, where :math:`f()`, :math:`g()` are rational polynomial function, and (P,L,H) normalized ground positions.
The rational function polynomial equation numerators and denominators each are 20-term cubic polynomial functions, which respects RPC00B convention.
Further details are given in `RPC in Geotiff`_ and `STDI-0002 2.1 (16Nov2000) specification document`_

Supported RPC format
--------------------

* DIMAP format (V1/V2)
* OSSIM keywordlist
* Geotiff RPC

Example
-------

The following parameter has to be set :
    * ``shareloc_home`` : Path to the Shareloc folder.

.. code-block:: bash

    import os
    from shareloc.rpc.rpc import RPC

    data_folder = os.path.join(shareloc_home,'valid')
    file_dimap = os.path.join(data_folder, f"rpc/RPC_PHR1B_P_201709281038393_SEN_PRG_FC_178609-001.XML")
    fctrat_dimap = RPC.from_any(file_dimap)


Direct location grids
=====================

Direct location grid is a sampled geometric models, which contains direct location at each grid cell for H0 to Hn altitude layers.
It can be viewed at 3D grid (row,col,h) as illustrated below :

.. figure:: images/direct_loc_multi_h.png
    :align: center
    :alt: direct location grid
    :width: 60%

    direct location grid

Shareloc grid format specifications
-----------------------------------

Shareloc grid must be a geotiff image, which contains 2 bands per altitude layer. One corresponding to x or longitude coordinates, the other corresponding to y, latitude coordinates

following metadata are needed

*  ALTITUDE_BX=Y : one per band X with altitude value Y
*  REF=EPSG:XXXX : coordinate reference system of ground coordinates

below an example of 9x5 grid composed of 3 altitude layers (-30m,485m,1000m). Each cell contains direct location at altitude layer of image position calculated from it's geotransform.

In the example below ``my_multi_h_grid`` is a 9x5x6 grid. ``my_multi_h_grid`` contains at index :math:`(row, col)` direct location
of :math:`((row + 0.5) * steprow + row0,  (col + 0.5) * stepcol + col0))`, for example with `(band, row, col)` convention
:math:`my\_multi\_h\_grid[0:1,1,2] = direct\_loc(row = 1250,col = 625,h = -30)`

.. code-block:: console

    $ gdalinfo my_multi_h_grid.tif

.. code-block:: console

    Driver: GTiff/GeoTIFF
    Files: test2.tif
    Size is 9, 5
    Coordinate System is `'
    Origin = (-312.500000000000000,-625.000000000000000)
    Pixel Size = (625.000000000000000,1250.000000000000000)
    Metadata:
      ALTITUDE_B0=-30.0
      ALTITUDE_B1=-30.0
      ALTITUDE_B2=485.0
      ALTITUDE_B3=485.0
      ALTITUDE_B4=1000.0
      ALTITUDE_B5=1000.0
      REF=EPSG:4326
    Image Structure Metadata:
      INTERLEAVE=PIXEL
    Corner Coordinates:
    Upper Left  (    -312.500,    -625.000)
    Lower Left  (    -312.500,    5625.000)
    Upper Right (    5312.500,    -625.000)
    Lower Right (    5312.500,    5625.000)
    Center      (    2500.000,    2500.000)
    Band 1 Block=9x5 Type=Float64, ColorInterp=Gray
    Band 2 Block=9x5 Type=Float64, ColorInterp=Undefined
    Band 3 Block=9x5 Type=Float64, ColorInterp=Undefined
    Band 4 Block=9x5 Type=Float64, ColorInterp=Undefined
    Band 5 Block=9x5 Type=Float64, ColorInterp=Undefined
    Band 6 Block=9x5 Type=Float64, ColorInterp=Undefined

Example
-------

The following parameter has to be set :
    * ``shareloc_home`` : Path to the Shareloc folder.

.. code-block:: bash

    import os
    from shareloc.grid import Grid
    data_folder = os.path.join(shareloc_home,'valid')
    eotiff_grid_path = os.path.join(data_folder, "ellipsoide", "loc_direct_grid_PHR_2013072139303958CP.tif")
    gri_geotiff = Grid(geotiff_grid_path)

References
__________

.. _`RPC in Geotiff`: http://geotiff.maptools.org/rpc_prop.html
.. _`STDI-0002 2.1 (16Nov2000) specification document`: http://geotiff.maptools.org/STDI-0002_v2.1.pdf