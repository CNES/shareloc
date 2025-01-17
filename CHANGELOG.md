# Changelog

Only the first "Unreleased" section of this file corresponding of next release can be updated along the development of each new, changed and fixed features.
When publication of a new release, the section "Unreleased" is blocked to the next chosen version and name of the milestone at a given date.
A new section Unreleased is opened then for next dev phase.

## Unreleased

### Added

### Changed

### Fixed

## 0.2.5 Margins for rectification grid (January 2025)

### Added

 - Add a configurable margin to rectification grid [#333]

## 0.2.4 Grid Interpolation method parameter (November 2024)

### Added

 - Interpolator type can be selected in rectification grid interpolation [#327]

### Fixed

- Fix remove interpolator to be cubic_legacy with scipy>=1.13 [#332] 
- Fix grid geomodel crash with Proj4 crs [#330]

## 0.2.3 Localisation grids for rectification (August 2024)

### Added

 - Rectification grids can be returned as localisation grids instead of displacement grids [#313]

### Changed

 - Allow nan values in geoid interpolation [#321]

### Fixed

 - Fix sphinx warnings during docs generation [#323]

## 0.2.2 M_PI identifier not found (July 2024)

### Fixed

- Fix M_PI identifier not found (Windows compatibility) [#322]

## 0.2.1 GRID model extrapolation and bug corrections, N views triangulation (June 2024)

### Added

- Grid extrapolation [#311]
- N views triangulation [#308]

### Fixed

- fix Numpy 2.0 pickle compatibility [#319]
- fix DTMIntersection failure with GRID models [#310]
- fix pylint failed [#315]
- fix DTMIntersection C++ binding interface [#316]

## 0.2.0 RPC, DTMintersection, Rectification C++ code optimisation (Mars 2024)

### Added

- Compute rectification grid by strips in C++  [#239, #247, #248, #249, #251, #270, #280, #297, #302]
- Force Y axis sign [#306]

### Changed

- Use RPC inverse coefficients by default in direct location [#287]
- Refactoring get_alt_offset() function [#296]

### Fixed

- fix altitude NaN filtering  [#210]
- fix image ROI outside the extent [#305]


## 0.2.0a2 Geoid interpolation correction, udpate on RPC C++ class (WIP), DTMintersection C++  (February 2024)

### Added

- RPCoptim C++ direct_loc_dtm()  method [#290]
- Epipolar angle rectification function [#294]
- DTMintersection C++ [#288, #295]

### Changed

- RPCoptim C++ geomodel class optimisation [#286]

### Fixed

- fix geoid interpolation  [#298]


## 0.2.0a1 GeoModel factory, RPC C++ class (WIP), by strip rectification, code optimisations (January 2024)

### Added

- Line of sight ending point [#240]
- Compute rectification grid by strips  [#236, #237, #238]
- GeoModel factory to instantiate geomodels [#221]
- Drive numba parallel by env SHARELOC_NUMBA_PARALLEL [#267]
- RPCoptim C++ geomodel class with direct loc h and inverse loc [#263, #225, #253, #276, #244, #254, #274, #273]
- Use rasterio rpc class to read image RPC [#23]

### Changed

- Code vectorization in DTMIntersection [#211]
- Use inverse geotransform as class argument to improve performances [#198]
- Remove useless new_axis [#226]
- Removed useless transforms in DTMIntersection [#215]
- RPC direct_loc_dtm refactoring [#275]

### Fixed

- fix build read the docs CI [#252]
- fix rpc geomodel_type typo in doc [#262]
- fix xarray pyarrow dependency [#291]

## 0.1.6 Minor bugs, dimap v3 experimental (October 2023)

### Added

- Add first dimap v3 experimental support [#155]

### Changed

- Clean logs for cars default output [#204]

### Fixed

- fix pandas.core.indexes.numeric evolution [#190]
- fix alt_min_max not referenced in los.py [#194]
- fix ValueError bug setting an array element with a sequence in py39 and py310 [#218]

## 0.1.5 Reliability enhancement (February 2023)

### Added

- Use scipy as ground truth for direct_loc_h [#120]
- Documentation in pre-commit [#169]
- Test with several python version through tox [#170]

### Changed

- Finalize grid files standard [#125]  
- Add epipolar footprint in prepare_rectification [#162]
- Remove netcdf4 package dependencies [#171, #180]

### Fixed

- Comment altitude extrapolation in RPC tests [#67]
- Analyze and comment margin in rectification and JP2 impact on Montpellier data [#86]
- Modified scipy version [#173]
- Clean code with black and pylint new version [#183] 

## 0.1.4 Clean packaging and documentation (December 2022)

### Added

- Authors file
- docstring sphinx autoapi generation in documentation

### Changed

- setup to python >=3.8

### Fixed

- Clean Makefile
- Clean python packaging
- fix math formula in sphinx generation
- fix recursion limit in astroid for sphinx autoapi and pylint
- fix lint errors (pylint, flake8, black, isort)

## 0.1.3 First Open Source Official Release - Quick fix 2 (April 2022)

### Fixed

- Documentation typos

## 0.1.2 First Open Source Official Release - Quick fix (March 2022)

### Fixed

- Documentation and quick start fixed (wrong URL and detailed install)

## 0.1.1 First Open Source Official Release (March 2022)

### Added

- Shareloc library first release
- geometric functions: localization, rectification, triangulation, earth elevation management
- geometric models: RPC and multi altitudes layers location grids
- Simple install with pip, clean dependencies (rasterio mainly) and Makefile for developper
- Documentation basics (README) and sphinx documentation to readthedocs
- Detailed pytest tests (tests directory)
- code quality (isort, black, flake8, pylint, sonarqube)
- Apache 2 license
