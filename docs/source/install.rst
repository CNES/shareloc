.. _install:

=======
Install
=======

Quick install
=============

* Create and activate virtualenv

.. code-block:: console

    $ python -m venv venv
    $ source venv/bin/activate
    $ python3 -m pip install --upgrade pip setuptools
    
* Install Shareloc from Pypi

.. code-block:: console

    $ pip install shareloc

Install from source
===================

* Clone Shareloc source code

.. code-block:: console

    $ git clone --depth 1 https://github.com/CNES/shareloc # For latest version
    $ git clone --depth 1 --branch LAST_TAG https://github.com/CNES/shareloc/ # For last stable version

* Install Shareloc from source

.. code-block:: console

    $ cd shareloc
    $ make install  # Shareloc is installed in `venv` directory

* Activate the virtual env

.. code-block:: console

    $ source venv/bin/activate