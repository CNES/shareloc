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


Numba parallelisation
=====================

**How to disable numba parallelisation ?**
By default Shareloc enables numba parallelisation. 
So if you want to work in single thread, set environment variable SHARELOC_NUMBA_PARALLEL to False.