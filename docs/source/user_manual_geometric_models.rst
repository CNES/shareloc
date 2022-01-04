.. _user_manual_geometric_models:


================
Geometric models
================

Shareloc handles two type of geometric models, RPC and Direct location grids. An upper layer class ``Localization`` can be used to handle both of them.

.. code-block:: console

    class Localization:
    """base class for localization function.
    Underlying model can be both multi layer localization grids or RPCs models
    """

    def __init__(self, model, elevation=None, image=None, epsg=None):
        """
        constructor
        :param model : geometric model
        :type model  : shareloc.grid or  shareloc.rpc
        :param elevation  : dtm or default elevation over ellipsoid if None elevation is set to 0
        :type elevation  : shareloc.dtm or float or np.ndarray
        :param image  : image class to handle geotransform
        :type image  : shareloc.image.image.Image
        :param epsg  : coordinate system of world points, if None model coordiante system will be used
        :type epsg  : int
        """

RPC
===

:term:`RPC` is analytics of ground (lon,lat,h) to image (r,c) mapping, it can be summarized as :math:`(r,c) = (f(P,L,H),g(P,L,H))`, where f(), g() are rational polynomial function, and (P,L,H) normalized ground positions.
The rational function polynomial equation numerators and denominators each are 20-term cubic polynomial functions, which respects RPC00B convention.
Further details are given in `RPC in Geotiff`_ and `STDI-0002 2.1 (16Nov2000) specification document`_

Supported RPC format
--------------------

* DIMAP format (V1/V2)
* OSSIM keywordlist
* Geotiff RPC

Example
-------

Direct location grids
=====================

Direct location grid is a sampled geometric models, which contains direct location at each grid cell for H0 .. Hn altitude layers.
It can be viewed at 3D grid (row,col,h) as illustrated below :


.. figure:: images/direct_loc_multi_h.png
    :align: center
    :alt: direct location grid
    :width: 60%

    direct location grid


Shareloc grid format
--------------------



Geotiff

metadata

GDAL infos

References
==========

.. _`RPC in Geotiff`: http://geotiff.maptools.org/rpc_prop.html
.. _`STDI-0002 2.1 (16Nov2000) specification document`: http://geotiff.maptools.org/STDI-0002_v2.1.pdf