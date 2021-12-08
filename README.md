

 
<div align="center">
  <a href="https://github.com/CNES/shareloc"><img src="docs/source/images/shareloc_picto.svg" alt="Shareloc" title="Shareloc"  width="20%"></a>

<h4>ShareLoc, a simple remote sensing geometric library</h4>

[![Python](https://img.shields.io/badge/python-v3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)](CONTRIBUTING.md)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0/)
[![Documentation](https://readthedocs.org/projects/shareloc/badge/?version=latest)](https://shareloc.readthedocs.io/?badge=latest)

<p>
  <a href="#overview">Overview</a> .
  <a href="#why-shareloc">Why Shareloc</a> .
  <a href="#quick-start">Quick Start</a> .
  <a href="#main-functions">Main Functions</a> .
  <a href="#documentation">Documentation</a> .
  <a href="#contribution">Contribution</a> .
</p>
</div>

## Overview

Shareloc is an open source remote sensing geolocation library. 

It performs image coordinates projections between sensor and ground and vice versa. 
Shareloc handles RPC models and direct location grids.


Direct localization at constant elevation |  Direct localization on DEM 
:-------------------------:|:-------------------------:
<img src="docs/source/images/shareloc_loc_ellipsoid.png" alt="drawing" width="400"/> |  <img src="docs/source/images/shareloc_loc_ellipsoid.png" alt="drawing" width="400"/>

Shareloc main functions : 

 * Direct/inverse localization at constant elevation (ellipsoidal earth model).
 * Direct localization on 2.5D DEM (w.r.t ellipsoid or geoid).
 * Line of sight triangulation.
 * Rectification grid creation.
 * Rectification grid interpolation.

## Why Shareloc

Shareloc development has been motivated by the need of a full python component for CNES studies, and the need of an underlying geometrical component for <a href="https://github.com/CNES/cars">CARS</a>.   


## Quick start

### Installation

Shareloc can be installed in a  [virtualenv](https://docs.python.org/3/library/venv) using the following commands:

```
git clone https://gitlab.cnes.fr/cars/shareloc.git
cd shareloc
make install
source venv/bin/activate # to go in installed dev environment
```

Dependencies : **git**, **make**

### Example

Shareloc is designed as a library. Please refer the [notebook directory](notebooks) for examples. 

## Documentation

Go to [Shareloc Main Documentation](https://shareloc.readthedocs.io/?badge=latest)  

## Contribution

To do a bug report or a contribution, see the [**Contribution Guide**](CONTRIBUTING.md).  
For project evolution, see [**Changelog**](CHANGELOG.md)