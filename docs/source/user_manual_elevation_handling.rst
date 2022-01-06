.. _user_manual_elevation_handling:


==================
Elevation handling
==================

The surface of the earth has to be modeled for geometrics functions like localization and rectification.
Two types of earth modelisation are available :

    * constant elevation over ellipsoid
    * :term:`DEM` which represents the surface of the earth

DEM
===

:term:`DEM` is a 2.5D representation of the surface of the earth.

Shareloc DEM constraints
------------------------

DEM must respects some constraints to be understandable by Shareloc :

 * format : DEM format has to be readable by GDAL (via ``rasterio``)
 * monolitic data : tiled DEM has to be mosaicked, using ``gdalbuildvrt`` command for example.
 * georeferenced : DEM must contains geotransform and :term:`CRS`.


.. code-block:: bash

    class DTM:

        def __init__(
            self,
            dtm_filename,
            geoid_filename=None,
            roi=None,
            roi_is_in_physical_space=True,
            fill_nodata=None,
            fill_value=0.0,
        ):




