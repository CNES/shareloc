.. _user_manual_functions:

=====================
Geolocation Functions
=====================

Shareloc has several main functions: 

- :ref:`Localization`
- :ref:`Triangulation`
- :ref:`Rectification Grid Computation`

Examples can be found in `tests directory <https://github.com/CNES/shareloc/tests/geofunctions>`_ of Shareloc source code.

Localization
============

Localization means mapping between image coordinates and ground ones using geometric model.

Localization functions can be called via its class or directly using geometric models. It can be simplified as the intersection of a :term:`LOS` and the terrain.

Following information are needed for localization functions:

 * **a geometric model**
 * **elevation information**, details in :ref:`user_manual_elevation_handling` section.
 * **image information** in order to use geotransform, details in :ref:`user_manual_conventions` section.
 * :term:`EPSG` **code** to specify ground coordinate system

``shareloc.geofunctions.localization.Localization`` class collect these data to set up the localization functions.
It is possible to use geometric model directly using ``shareloc.grid.Grid`` and ``shareloc.rpc.rpc.RPC`` for advanced used.

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
        :type image  : shareloc.image.Image
        :param epsg  : coordinate system of world points, if None model coordinate system will be used
        :type epsg  : int
        """


Direct Localization
-------------------

Direct localization returns ground coordinates  :math:`(\lambda,\phi,h)` for image position (row, column) :math:`(\lambda,\phi) = direct\_localization(row,column)`.

:math:`h` can be explicitly set as input (in case of constant altitude over ellipsoid) : :math:`(\lambda,\phi,h) = direct\_localization(row,column,h)`, or found via the underlying :term:`DEM` : :math:`(\lambda,\phi) = direct\_localization(row,column)`

:math:`(\lambda,\phi)` can be geographic (lat,lon) or cartographic depending on the chosen :term:`CRS`

.. code-block:: bash

    def direct(self, row, col, h=None, using_geotransform=False):
        """
        direct localization
        :param row :  sensor row
        :type row : float
        :param col : sensor col
        :type col : float
        :param h: altitude, if none DTM is used
        :type h : float
        :param using_geotransform: using_geotransform
        :type using_geotransform : boolean
        :return coordinates : [lon,lat,h] (3D np.array)
        """


Inverse Localization
--------------------

inverse localization returns image position (row,column) for ground coordinates :math:`(\lambda,\phi,h)`  :math:`(row,col) = inverse\_localization(\lambda,\phi,h)`.

.. code-block:: bash

    def inverse(self, lon, lat, h=None, using_geotransform=False):
        """
        inverse localization
        :param lat :  latitude (or y)
        :param lon : longitude (or x)
        :param h : altitude
        :param using_geotransform: using_geotransform
        :type using_geotransform : boolean
        :return coordinates : [row,col,h] (2D np.array)
        :rtype numpy.array
        """


Colocalization
--------------

colocalization returns image positions (row2,col2) in image 2 from (row1,col1) position in image 1

.. code-block:: bash

    def coloc(model1, model2, row, col, elevation=None, image1=None, image2=None, using_geotransform=False):
        """
        Colocalization : direct localization with model1, then inverse localization with model2

        :param model1: geometric model 1
        :type model1: shareloc.grid or  shareloc.rpc
        :param model2: geometric model 2
        :type model2: shareloc.grid or  shareloc.rpc
        :param row: sensor row
        :type row: int or 1D numpy array
        :param col: sensor col
        :type col: int or 1D numpy array
        :param elevation: elevation
        :type elevation: shareloc.dtm or float or 1D numpy array
        :param image1  : image class to handle geotransform
        :type image1  : shareloc.image.Image
        :param image2  : image class to handle geotransform
        :type image2  : shareloc.image.Image
        :param using_geotransform: using_geotransform
        :type using_geotransform : boolean
        :return: Corresponding sensor position [row, col, True] in the geometric model 2
        :rtype : Tuple(1D np.array row position, 1D np.array col position, 1D np.array True)
        """


Triangulation
=============

Triangulation gives 3D intersections between :term:`LOS` coming from 2 geometric models.

Triangulation is calculated according to the following formula:

:math:`x= \left(\sum_i I-\hat v_i \hat v_i^\top\right)^{-1} \left(\sum_i (I-\hat v_i \hat v_i^\top) s_i\right)`

where :math:`v_i` is the orientation of the :term:`LOS` i and :math:`s_i` the hat of the :term:`LOS` i

.. code-block:: bash

    def sensor_triangulation(
        matches,
        geometrical_model_left,
        geometrical_model_right,
        left_min_max=None,
        right_min_max=None,
        residues=False,
        fill_nan=False,
    ):
        """
        triangulation in sensor geometry

        according to the formula:
        .. math::
            x =
            \\left(\\sum_i I-\\hat v_i \\hat v_i^\\top\\right)^{-1} \\left(\\sum_i (I-\\hat v_i \\hat v_i^\\top) s_i\\right)
        Delvit J.M. et al. "The geometric supersite of Salon de Provence", ISPRS Congress Paris, 2006.


        :param matches :  matches in sensor coordinates Nx[row (left), col (left), row (right), col (right)]
        :type matches : np.array
        :param geometrical_model_left : left image geometrical model
        :type geometrical_model_left : shareloc.grid or shareloc.rpc
        :param geometrical_model_right : right image geometrical model
        :type geometrical_model_right : shareloc.grid or shareloc.rpc
        :param left_min_max : left min/max for los creation, if None model min/max will be used
        :type left_min_max : list
        :param right_min_max : right min/max for los creation, if None model min/max will be used
        :type right_min_max : list
        :param residues : calculates residues (distance in meters between los and 3D points)
        :type residues : boolean
        :param fill_nan : fill numpy.nan values with lon and lat offset if true (same as OTB/OSSIM), nan is returned
            otherwise
        :type fill_nan : boolean
        :return intersections in cartesian crs, intersections in wgs84 crs and optionnaly residues
        :rtype (numpy.array,numpy,array,numpy.array)
        """

References :
------------

- Delvit J.M. et al. **The geometric supersite of Salon de Provence**, ISPRS Congress Paris, 2006. (`http://isprs.free.fr/documents/Papers/T11-50.pdf <http://isprs.free.fr/documents/Papers/T11-50.pdf>`_)


Rectification Grid Computation
==============================

:term:`Rectification` or stereo-rectification refers to the image transformation in epipolar geometry.

A rectification grid is a displacement grid used to resample sensor gemetry to epipolar one.
Shareloc rectification grids respects OTB convention for displacement grids. 

To generate the images in epipolar geometry from the grids computed by shareloc and the original images, one can refer to the Orfeo Toolbox documentation `here <https://www.orfeo-toolbox.org/CookBook/recipes/stereo.html#resample-images-in-epipolar-geometry>`_ .
Algorithm details can be found in reference below.

.. code-block:: bash

    def compute_stereorectification_epipolar_grids(
        left_im, geom_model_left, right_im, geom_model_right, elevation=0.0, epi_step=1, elevation_offset=50.0
    ):
        """
        Compute stereo-rectification epipolar grids

        :param left_im: left image
        :type left_im: shareloc.image object
        :param geom_model_left: geometric model of the left image
        :type geom_model_left: shareloc.grid or  shareloc.rpc
        :param right_im: right image
        :type right_im: shareloc.image object
        :param geom_model_right: geometric model of the right image
        :type geom_model_right: shareloc.grid or  shareloc.rpc
        :param elevation: elevation
        :type elevation: shareloc.dtm or float
        :param epi_step: epipolar step
        :type epi_step: int
        :param elevation_offset: elevation difference used to estimate the local tangent
        :type elevation_offset: float
        :return: return :
            - left epipolar grid, shareloc.image object convention [[row displacement, col displacement], nb rows, nb cols]
            - right epipolar grid, shareloc.image object convention [[row displacement, col displacement], nb rows, nb cols]
            - number of rows of the epipolar image, int
            - number of columns of the epipolar image, int
            - mean value of the baseline to sensor altitude ratio, float
        :rtype: Tuple
        """


References :
------------
- Youssefi D., Michel, J., Sarrazin, E., Buffe, F., Cournet, M., Delvit, J.,  L'Helguen, C., Melet, O., Emilien, A., Bosman, J., 2020. **CARS: A photogrammetry pipeline using dask graphs to construct a global 3d model**. IGARSS - IEEE International Geoscience and Remote Sensing Symposium.(`https://ieeexplore.ieee.org/document/9324020 <https://ieeexplore.ieee.org/document/9324020>`_)