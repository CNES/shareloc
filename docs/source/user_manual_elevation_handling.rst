.. _user_manual_elevation_handling:


==================
Elevation handling
==================

The surface of the earth has to be modeled for geometrics functions localization and rectification.
Two types of earth model are available :

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

:term:`DEM` can be used for Localization on DEM function and Rectification using `shareloc.dtm.DTM` class.

Since Shareloc works w.r.t elllipsoid by default, geoid height has to be removed from :term:`DEM` if w.r.t geoid.
Shareloc provides egm96_15 in its  `data <https://gitlab.cnes.fr/cars/shareloc/-/blob/master/valid/dtm/geoid/egm96_15.gtx>`

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

For example, the `SRTM <https://www2.jpl.nasa.gov/srtm/>`_ data corresponding to the zone to process can be used through the `otbcli_DownloadSRTMTiles <https://www.orfeo-toolbox.org/CookBook/Applications/app_DownloadSRTMTiles.html>`_ OTB command.
