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
The trasnformation between index and physical point use the geotransform of the iamge.

``shareloc.image.image.Image`` class manage the conversion between index and physical point.


It allow to manage easily images :term:`ROI`.

(see example)
(ref OTB software guide chapter 5 p62)



