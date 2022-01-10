.. _install:

=======
Install
=======
.. _dependencies:

Dependencies
=============

Shareloc depends on few python packages

Quick install
=============

* Clone Shareloc Source code

.. code-block:: console

    $ git clone --depth 1 https://gitlab.cnes.fr/cars/shareloc.git # For latest version
    $ git clone --depth 1 --branch LAST_TAG https://gitlab.cnes.fr/cars/shareloc.git # For last stable version

* Install Shareloc

.. code-block:: console

    $ cd shareloc
    $ make install  # Shareloc is installed in `venv` directory

* Activate the virtual env

.. code-block:: console

    $ source venv/bin/activate

Advanced Install
================
The following steps are defined in `Makefile`  ``install`` command.

Virtualenv
----------
First create a virtualenv and upgrade main pip packages.

.. code-block:: console

    $ virtualenv -p python venv/
    $ source venv/bin/activate
    $ python3 -m pip install --upgrade pip setuptools

Required python packages
------------------------

Shareloc python package requires some python packages to be installed before:

* **numpy**, **cython**: They have to be installed separately otherwise some dependencies won't be correctly installed.
* **rasterio**: On some systems, it has to be installed from source to fit local GDAL version.

Here are the correspondent commands to install these prior dependencies:

.. code-block:: console

    $ virtualenv -p python venv/
    $ source venv/bin/activate
    $ python3 -m pip install --upgrade cython numpy
    $ python3 -m pip install --no-binary rasterio rasterio
    $ python3 -m pip install pygdal=="$(gdal-config --version).*"


Environment variables
---------------------

TEST_PATH must be set at shareloc root dir (TEST_PATH is going to be suppressed)