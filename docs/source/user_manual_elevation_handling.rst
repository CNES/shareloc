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

:term:`DEM` can be used for Localization on DEM function and Rectification using `shareloc.geofunctions.DTMIntersection` class.

Since Shareloc works w.r.t elllipsoid by default, geoid height has to be removed from :term:`DEM` if w.r.t geoid.

Shareloc provides egm96_15 in its  `data <https://raw.githubusercontent.com/CNES/shareloc/tests/data/dtm/geoid/egm96_15.gtx>`_

.. code-block:: bash

    class DTMdtm_reader(Image):

        def __init__(
        	self,
        	dtm_filename,
        	geoid_filename=None,
        	read_data=False,
        	roi=None,
        	roi_is_in_physical_space=False,
        	fill_nodata="rio_fillnodata",
        	fill_value=None,
        ):

For example, the `SRTM <https://www2.jpl.nasa.gov/srtm/>`_ data corresponding to the zone to process can be used through the `otbcli_DownloadSRTMTiles <https://www.orfeo-toolbox.org/CookBook/Applications/app_DownloadSRTMTiles.html>`_ OTB command.

Limitations
-----------

 * nodata are not (yet) handled in `shareloc.geofunctions.DTMIntersection` code. Thus a filling strategy should be set when using nodata DTM. This can be done by setting `fill_nodata` arg in `shareloc.geofunctions.DTMIntersection`. Filling strategy examples can be found in test `test_dtm_image`

