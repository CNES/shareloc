.. _faq:

===============
FAQ
===============

Input Pleiades data
===================

**Why does my rectification not work on the Pleiades product ?**

If you use a Pleiades product, make sure to use dimap as an input and not JP2.
Because the size of the pixels is negative in the Y dimension.
This note is available to all tools like CARS or OTB.

RPC formats
===========

**How can i transform a RPC .geom file into a gdal geotiff with RPC tags ?**

a tool ``shareloc-rpcconverter`` is available to convert rpc into geotiff or geotiff with RPB.
The following command is an example to add `rpc.geom` rpc coefficient to `img.tif`

.. code-block:: console

    shareloc-rpcconverter rpc.geom img.tif geotiff

.. code-block:: console

    usage: shareloc-rpcconverter [-h] [--override] [--rpc_format RPC_FORMAT] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] input_rpc output_file

    convert RPC models

    positional arguments:
      input_rpc             input rpc
      output_file           output_file

    options:
      -h, --help            show this help message and exit
      --override            override geotiff rpc if already present.
      --rpc_format RPC_FORMAT
                            output rpc format in geotiff, geotiff_rpb, rpb, json. if format is geotiff or geotiff_rpb output_file must be an existing tif, which will be
                            updated with rpc tifftag,or external .RPB file in case of geotiff_rpb, if rpb, output_file must be a .RPB.
      --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Logger level (default: INFO. Should be one of (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    This script takes as input an RPC model and transforms it to another format"




Numba parallelisation
=====================

**How to disable numba parallelisation ?**

By default Shareloc enables numba parallelisation. 
So if you want to work in single thread, set environment variable SHARELOC_NUMBA_PARALLEL to False.


