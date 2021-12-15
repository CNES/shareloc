.. _glossary:

========
Glossary
========

Shareloc common words and shortened terms are detailed here.

To update, follow `glossary sphinx documentation`_ in RST source documentation.

.. glossary::
    CARS
      means CNES Algorithms to Reconstruct Surface (ou Chaîne Automatique de Restitution Stéréoscopique en français)

    DEM
      `Digital Elevation Model`_. Usually means all elevation models in raster: DSM, DTM,...

    DSM
      Digital Surface Model. Represents the earth's surface and includes all objects on it.
      CARS generates DSMs. See `Digital Elevation Model`_

    DTM
      Digital Terrain Model. Represents bare ground surface without any objects like plants and buildings
      You need another tool to generate DTM from CARS DSM. See `Digital Elevation Model`_

    rectification
      `Image rectification`_ is a transformation process used to project images onto a common image plane.
      In CARS, the epipipolar geometry rectification is used.

    ROI
      `Region of Interest`_ means a subpart of the `DSM` raster in CARS.
      It can be defined by a file or a bounding box.

.. _`Digital Elevation Model`: https://en.wikipedia.org/wiki/Digital_elevation_model
.. _`Digital Surface Model`: https://en.wikipedia.org/wiki/Digital_elevation_model
.. _`epipolar geometry`: https://en.wikipedia.org/wiki/Epipolar_geometry
.. _`Image rectification`: https://en.wikipedia.org/wiki/Image_rectification
.. _`Region of Interest`: https://en.wikipedia.org/wiki/Region_of_interest

.. _`glossary sphinx documentation`: https://sublime-and-sphinx-guide.readthedocs.io/en/latest/glossary.html
