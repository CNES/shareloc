# Changelog

Only the first "Unreleased" section of this file corresponding of next release can be updated along the development of each new, changed and fixed features.
When publication of a new release, the section "Unreleased" is blocked to the next chosen version and name of the milestone at a given date.
A new section Unreleased is opened then for next dev phase.

## Unreleased

### Added

### Changed

### Fixed


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
