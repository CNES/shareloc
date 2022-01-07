.. _user_manual_conventions:


This sections details Shareloc convention for images and geometric models.

================
Pixel convention
================

Geometric model convention for coordinates is, by default, [0.5,0.5] at the center of the first pixel of the image, as illustrated below.

.. figure:: images/convention_pixel.png
    :align: center
    :alt: convention pixel
    :width: 60%

    convention pixel

When dealing with ``shareloc.rpc.rpc.RPC`` it is possible to fix the center at [0,0] by setting the ``topleftconvention`` option to ``False``.

.. code-block:: bash

    @classmethod
    def from_any(cls, primary_file, secondary_file=None, topleftconvention=True):


================
Image convention
================

When dealing with image we use ``index`` to refer the position in the 2D array and ``Physical point`` the position of the pixel.
The transformation between index and physical point use the geotransform of the image. `OTB Software guide chapter 5.1`_ gives complete descriptions of these concepts.

``shareloc.image.image.Image`` class manage the conversion between index and physical point.

It allow to manage easily images :term:`ROI` without any processing on geomtric model. See `CARS Faq <https://cars.readthedocs.io/en/latest/faq.html#faq>` or `GDAL translate command <https://gdal.org/programs/gdal_translate.html>` to create image extract.

see `image example <https://github.com/CNES/shareloc/tests/test_image.py>`



.. _`OTB Software guide chapter 5.1` : https://www.orfeo-toolbox.org/packages/archives/Doc/SoftwareGuide-6.6.0.pdf
