# see: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
# Dependencies : python3 venv
# Some Makefile global variables can be set in make command line
# Recall: .PHONY  defines special targets not associated with files

############### GLOBAL VARIABLES ######################
.DEFAULT_GOAL := help
# Set shell to BASH
SHELL := /bin/bash

# Set LOGLEVEL if not defined in command line
# Example: LOGLEVEL="DEBUG" make help
ifndef LOGLEVEL
	LOGLEVEL = "INFO"
endif

# Set Virtualenv directory name
# Exemple: VENV="other-venv/" make install
ifndef VENV
	VENV = "venv"
endif

# Check shareloc install
CHECK_SHARELOC = $(shell ${VENV}/bin/python -m pip list|grep shareloc)


# Get setuptools_scm version 
VERSION = $(shell python3 -c 'from shareloc import __version__; print(__version__)')
VERSION_MIN = $(shell echo ${VERSION} | cut -d . -f 1,2,3)

# Browser definition
define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT
BROWSER := ${VENV}/bin/python -c "$$BROWSER_PYSCRIPT"

# Python global variables definition
PYTHON_VERSION_MIN = 3.9


# Set PYTHON if not defined in command line
# Example: PYTHON="python3.10" make venv to use python 3.10 for the venv
# By default the default python3 of the system.
ifndef PYTHON
	PYTHON = "python3"
endif
PYTHON_CMD=$(shell command -v $(PYTHON))

PYTHON_VERSION_CUR=$(shell $(PYTHON_CMD) -c 'import sys; print("%d.%d"% sys.version_info[0:2])')
PYTHON_VERSION_OK=$(shell $(PYTHON_CMD) -c 'import sys; cur_ver = sys.version_info[0:2]; min_ver = tuple(map(int, "$(PYTHON_VERSION_MIN)".split("."))); print(int(cur_ver >= min_ver))')

############### Check python version supported ############

ifeq (, $(PYTHON_CMD))
    $(error "$(PYTHON_CMD) not found in $(PATH)")
endif

ifeq ($(PYTHON_VERSION_OK), 0)
    $(error "Requires python version >= $(PYTHON_VERSION_MIN). Current version is $(PYTHON_VERSION_CUR)")
endif

################ MAKE targets by sections ######################

.PHONY: help
help: ## this help
	@echo "      SHARELOC MAKE HELP"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' | sort

## Install section

.PHONY: venv
venv: ## create virtualenv in "venv" dir if not exists
	@test -d ${VENV} || $(PYTHON_CMD) -m venv ${VENV}
	@${VENV}/bin/python -m pip install --upgrade pip setuptools wheel # no check to upgrade each time install is run
	@touch ${VENV}/bin/activate

.PHONY: install
install: venv  ## install the package in dev mode in virtualenv
	@[ "${CHECK_SHARELOC}" ] || echo "Install shareloc package from local directory"
	@[ "${CHECK_SHARELOC}" ] || ${VENV}/bin/python -m pip install -e .[dev,docs,notebook]
	@test -f .git/hooks/pre-commit || echo "Install pre-commit"
	@test -f .git/hooks/pre-commit || ${VENV}/bin/pre-commit install -t pre-commit
	@test -f .git/hooks/pre-push || ${VENV}/bin/pre-commit install -t pre-push
	@echo "Shareloc ${VERSION} installed in dev mode in virtualenv ${VENV} with documentation"
	@echo " virtualenv usage : source ${VENV}/bin/activate; python3 -c 'import shareloc'"

## Test section

.PHONY: test
test: ## run tests and coverage quickly with the default Python
	@${VENV}/bin/pytest -o log_cli_level=${LOGLEVEL} --cov-config=.coveragerc --cov --cov-report=term-missing

.PHONY: test-all
test-all: ## run tests on every Python version with tox
	@${VENV}/bin/tox -r -p auto  ## recreate venv (-r) and parallel mode (-p auto)

.PHONY: coverage
coverage: ## check code coverage quickly with the default Python
	@${VENV}/bin/coverage run --source shareloc -m pytest
	@${VENV}/bin/coverage report -m
	@${VENV}/bin/coverage html
	$(BROWSER) htmlcov/index.html

## Code quality, linting section

### Format with isort and black

.PHONY: format
format: format/isort format/black  ## run black and isort formatting

.PHONY: format/isort
format/isort: ## run isort formatting
	@echo "+ $@"
	@${VENV}/bin/isort shareloc tests

.PHONY: format/black
format/black: ## run black formatting
	@echo "+ $@"
	@${VENV}/bin/black shareloc tests

### Check code quality and linting : isort, black, flake8, pylint

.PHONY: lint
lint:  lint/isort lint/black lint/flake8 lint/pylint ## check code quality and linting

.PHONY: lint/isort
lint/isort: ## check imports style with isort
	@echo "+ $@"
	@${VENV}/bin/isort --check shareloc tests

.PHONY: lint/black
lint/black: ## check global style with black
	@echo "+ $@"
	@${VENV}/bin/black --check shareloc tests

.PHONY: lint/flake8
lint/flake8: ## check linting with flake8
	@echo "+ $@"
	@${VENV}/bin/flake8 shareloc tests

.PHONY: lint/pylint
lint/pylint: ## check linting with pylint
	@echo "+ $@"
	@set -o pipefail; ${VENV}/bin/pylint shareloc tests --rcfile=.pylintrc --output-format=parseable | tee pylint-report.txt # pipefail to propagate pylint exit code in bash

.PHONY: lint/mypy
lint/mypy: ## check linting type hints with mypy
	@echo "+ $@"
	@${VENV}/bin/mypy shareloc tests

## Documentation section

.PHONY: docs
docs: install ## generate Sphinx HTML documentation, including API docs
	@${VENV}/bin/sphinx-build -M clean docs/source/ docs/build
	@${VENV}/bin/sphinx-build -M html docs/source/ docs/build -W --keep-going
	$(BROWSER) docs/build/html/index.html

## Notebook section

.PHONY: notebook
notebook: ## Install Jupyter notebook kernel with venv
	@echo "Install Jupyter Kernel and launch Jupyter notebooks environment"
	@${VENV}/bin/python -m ipykernel install --sys-prefix --name=shareloc-$(VENV) --display-name=shareloc-$(VERSION)
	@echo " --> After virtualenv activation, please use following command to launch local jupyter notebook to open Notebooks:"
	@echo "jupyter notebook"

.PHONY: notebook-clean-output ## Clean Jupyter notebooks outputs
notebook-clean-output:
	@echo "Clean Jupyter notebooks"
	jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace notebooks/*.ipynb

## Docker section

docker: ## Build docker image (and check Dockerfile)
	@echo "Check Dockerfile with hadolint"
	@docker pull hadolint/hadolint
	@docker run --rm -i hadolint/hadolint < Dockerfile
	@echo "Build Docker image shareloc ${VERSION_MIN}"
	@docker build -t cnes/shareloc:${VERSION_MIN} -t cnes/shareloc:latest .

## Release section

.PHONY: dist
dist: clean install ## clean, install, builds source and wheel package
	@${VENV}/bin/python -m pip install --upgrade build
	@${VENV}/bin/python -m build
	ls -l dist

.PHONY: release
release: dist ## package and upload a release
	@twine check dist/*
	@twine upload dist/* --verbose ##  update your .pypirc accordingly

## Clean section

.PHONY: clean
clean: clean-venv clean-build clean-precommit clean-pyc clean-test clean-docs clean-notebook ## remove all venv, build, precommit, python artifacts, test, docs, notebook cache

.PHONY: clean-venv
clean-venv: ## clean installed venv
	@echo "+ $@"
	@rm -rf ${VENV}

.PHONY: clean-build
clean-build: ## clean build artifacts
	@echo "+ $@"
	@rm -fr build/
	@rm -fr dist/
	@rm -fr .eggs/
	@find . -name '*.egg-info' -exec rm -fr {} +
	@find . -name '*.egg' -exec rm -f {} +

.PHONY: clean-precommit
clean-precommit: ## clean precommit
	@rm -f .git/hooks/pre-commit
	@rm -f .git/hooks/pre-push

.PHONY: clean-pyc
clean-pyc: ## clean python file artifacts
	@echo "+ $@"
	@find . -type f -name "*.py[co]" -exec rm -fr {} +
	@find . -type d -name "__pycache__" -exec rm -fr {} +
	@find . -name '*~' -exec rm -fr {} +

.PHONY: clean-test
clean-test: ## clean test and coverage artifacts
	@echo "+ $@"
	@rm -fr .tox/
	@rm -f .coverage
	@rm -rf .coverage.*
	@rm -rf coverage.xml
	@rm -fr htmlcov/
	@rm -fr .pytest_cache
	@rm -f pytest-report.xml
	@rm -f pylint-report.txt
	@rm -f debug.log

.PHONY: clean-docs
clean-docs: ## clean generated docs
	@echo "+ $@"
	@rm -rf docs/build/
	@rm -rf docs/source/api_reference/
	@rm -rf docs/source/apidoc/

.PHONY: clean-notebook
clean-notebook: ## clean notebooks cache
	@echo "+ $@"
	@find . -type d -name ".ipynb_checkpoints" -exec rm -fr {} +

.PHONY: clean-docker
clean-docker: ## clean shareloc docker images
		@echo "+ $@"
		@echo "Clean Docker image shareloc ${VERSION_MIN}"
		@docker image rm cnes/shareloc:${VERSION_MIN}
		@docker image rm cnes/shareloc:latest
