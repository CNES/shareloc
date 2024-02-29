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
 * EPSG:4326 only if you use DTMIntersection(c++) with RPCoptim
 
Since Shareloc works w.r.t elllipsoid by default, geoid height has to be removed from :term:`DEM` if w.r.t geoid.

Shareloc provides `egm96_15` in its [data](https://github.com/CNES/shareloc/tree/master/tests/data/dtm/geoid).

The term :term:`DEM` can be used for Localization on DEM function and Rectification using the `shareloc.geofunctions.DTMIntersection` class. This class is initialized with a `dtm_reader` object as follows:

.. code-block:: Python

    import bindings_cpp
    
    dtm_image = dtm_reader( # dtm_reader with default arguments
        dtm_filename,
        geoid_filename=None,
        roi=None,
        roi_is_in_physical_space=False,
        fill_nodata="rio_fillnodata",
        fill_value=None,
    )
    dtm_py = DTMIntersection(#python version of DTMIntersection
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )
    dtm_cpp = bindings_cpp.DTMIntersection(#C++ version of DTMIntersection
        dtm_image.epsg,
        dtm_image.alt_data,
        dtm_image.nb_rows,
        dtm_image.nb_columns,
        dtm_image.transform,
    )



For example, the `SRTM <https://www2.jpl.nasa.gov/srtm/>`_ data corresponding to the zone to process can be used through the `otbcli_DownloadSRTMTiles <https://www.orfeo-toolbox.org/CookBook/Applications/app_DownloadSRTMTiles.html>`_ OTB command.

Limitations
-----------

 * nodata are not (yet) handled in `shareloc.geofunctions.DTMIntersection` code. Thus a filling strategy should be set when using nodata DTM. This can be done by setting `fill_nodata` arg in `shareloc.geofunctions.DTMIntersection`. Filling strategy examples can be found in test `test_dtm_image`


