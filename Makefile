# Autodocumented Makefile
# see: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html

## INPUT VARIABLES ##

# Set LOGLEVEL if not defined in command line
# Example: LOGLEVEL="DEBUG" make help
ifndef LOGLEVEL
	LOGLEVEL = "INFO"
endif

## GLOBAL VARIABLES ##

# Set shell to BASH
SHELL := /bin/bash

# Set Virtualenv directory name
ifndef VENV
	VENV = "venv"
endif

# SHARELOC variables
CHECK_SHARELOC = $(shell ${VENV}/bin/python -m pip list|grep shareloc)
SHARELOC_VERSION = $(shell python3 setup.py --version)

# MAKE TARGETS DEFINITION
.PHONY: help check venv install test test-ci lint format notebook clean

help: ## this help
	@echo "      SHARELOC MAKE HELP  LOGLEVEL=${LOGLEVEL}"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

venv: ## create virtualenv in "venv" dir if not exists
	@test -d ${VENV} || python3 -m venv ${VENV}
	@${VENV}/bin/pip install --upgrade pip setuptools # no check to upgrade each time
	@touch ${VENV}/bin/activate

install: venv ## install shareloc in dev mode
	@[ "${CHECK_SHARELOC}" ] || ${VENV}/bin/pip install --verbose -e .[dev]
	@test -f .git/hooks/pre-commit || echo "  Install pre-commit hook"
	@test -f .git/hooks/pre-commit || ${VENV}/bin/pre-commit install -t pre-commit
	@echo "Shareloc ${SHARELOC_VERSION} installed in dev mode in virtualenv ${VENV}"
	@echo "Shareloc venv usage : source ${VENV}/bin/activate; python3 -c 'import shareloc'"

install-doc: install  ## install shareloc with Sphinx documentation dependencies
	@[ "${CHECK_SHARELOC}" ] || ${VENV}/bin/pip install --verbose -e .[doc]
	@echo "Shareloc ${SHARELOC_VERSION} in virtualenv ${ENV} installed with Sphinx docs dependencies"

doc: install-doc ## build sphinx documentation
	@${VENV}/bin/sphinx-build -M clean docs/source/ docs/build
	# @${VENV}/bin/sphinx-apidoc -o docs/source/apidoc/ shareloc
	@${VENV}/bin/sphinx-build -M html docs/source/ docs/build


test: install ## run all tests + coverage html
	@${VENV}/bin/pytest -o log_cli=true -o log_cli_level=${LOGLEVEL} --cov-config=.coveragerc --cov-report html --cov

test-ci: install ## run all + coverage for ci to sonarqube
	@${VENV}/bin/pytest -o log_cli=true -o log_cli_level=${LOGLEVEL} --junitxml=pytest-report.xml --cov-config=.coveragerc --cov-report xml --cov

lint: install ## run lint tools (depends install)
	@echo "Linting isort check"
	@${VENV}/bin/isort --check shareloc tests
	@echo "Linting black check"
	@${VENV}/bin/black --check shareloc tests
	@echo "Linting flake8 check"
	@${VENV}/bin/flake8 shareloc tests
	@echo "Linting pylint check"
	@set -o pipefail; ${VENV}/bin/pylint shareloc tests --rcfile=.pylintrc --output-format=parseable | tee pylint-report.txt # pipefail to propagate pylint exit code in bash

format: install  ## run black and isort formatting (depends install)
	@${VENV}/bin/isort shareloc tests
	@${VENV}/bin/black shareloc tests

notebook: install ## Install Jupyter notebook kernel with venv and shareloc install
	@echo "Install Jupyter Kernel and launch Jupyter notebooks environment"
	@${VENV}/bin/python -m ipykernel install --sys-prefix --name=shareloc-$(VENV) --display-name=shareloc-$(SHARELOC_VERSION)
	@echo "--> After Shareloc virtualenv activation, please use following command to launch local jupyter notebook to open Shareloc Notebooks:"
	@echo "jupyter notebook"

clean: ## clean: remove all generated files: venv, cache, ...
	@rm -rf ${VENV}
	@rm -rf dist
	@rm -rf shareloc.egg-info
	@find . -type d -name __pycache__ -exec rm -r {} \+
	@rm -rf .eggs
	@rm -f .coverage
	@rm -rf .coverage.*
