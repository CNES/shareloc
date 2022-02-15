========
Glossary
========

Shareloc common words and shortened terms are detailed here.

.. glossary::
    CARS
      means CNES Algorithms to Reconstruct Surface (ou Chaîne Automatique de Restitution Stéréoscopique en français)
      See `CARS Documentation`_

    CRS
      `Coordinate Reference System`_ characterize location on earth.


    DEM
      `Digital Elevation Model`_. Usually means all elevation models in raster: DSM, DTM,...

    DSM
      Digital Surface Model. Represents the earth's surface and includes all objects on it.
      CARS generates DSMs. See `Digital Elevation Model`_
      
    DTM
      Digital Terrain Model. Represents the earth's surface without the objects. 
    
    EPSG
      EPSG Geodetic Parameter Dataset (also EPSG registry) is a public registry of geodetic datums, spatial reference systems, Earth ellipsoids, coordinate transformations and related units of measurement. 
      See `EPSG`_
      
    LOS
      Line of Sight is the viewing ray at a given position of the sensor.

    Rectification
      `Image rectification`_ is a transformation process used to project images onto a common image plane.
      In CARS, the epipipolar geometry rectification is used.

    RPC
      `Rational polynomial coefficient`_ is analytics compact representation of ground to image mapping.

    ROI
      `Region of Interest`_ means a subpart of the data.
      It can be defined by it's extent.

.. _`Digital Elevation Model`: https://en.wikipedia.org/wiki/Digital_elevation_model
.. _`Digital Surface Model`: https://en.wikipedia.org/wiki/Digital_elevation_model
.. _`epipolar geometry`: https://en.wikipedia.org/wiki/Epipolar_geometry
.. _`Image rectification`: https://en.wikipedia.org/wiki/Image_rectification
.. _`Region of Interest`: https://en.wikipedia.org/wiki/Region_of_interest
.. _`Rational polynomial coefficient`: https://en.wikipedia.org/wiki/Rational_polynomial_coefficient
.. _`glossary sphinx documentation`: https://sublime-and-sphinx-guide.readthedocs.io/en/latest/glossary.html
.. _`Coordinate Reference System`: https://en.wikipedia.org/wiki/Spatial_Reference_system
.. _`CARS Documentation`: https://cars.readthedocs.io/
.. _`EPSG`: https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset
