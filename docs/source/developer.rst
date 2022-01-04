================
Developer Manual
================

Shareloc is an open source software : don't hesitate to hack it and contribute !

Please go to `the GitHub repository`_  for source code.

Read also `Shareloc Contribution guide`_ with `LICENCE <https://raw.githubusercontent.com/CNES/shareloc/master/LICENSE>`_ and `Contributor Licence Agrements <https://github.com/CNES/shareloc/tree/master/docs/source/CLA>`_.

**Contact:** cars AT cnes.fr

Developer Install
=================
:ref:`Install` procedure is followed.

The detailed development install method is described in `Makefile`

Particularly, it uses the following pip editable install:

.. code-block:: console

    pip install -e .[dev]

With this pip install mode, source code modifications directly impacts Shareloc code.

Coding guide
============

Here are some rules to apply when developing a new functionality:

* **Comments:** Include a comments ratio high enough and use explicit variables names. A comment by code block of several lines is necessary to explain a new functionality.
* **Test**: Each new functionality shall have a corresponding test in its module's test file. This test shall, if possible, check the function's outputs and the corresponding degraded cases.
* **Documentation**: All functions shall be documented (object, parameters, return values).
* **Use type hints**: Use the type hints provided by the `typing` python module.
* **Use doctype**: Follow sphinx default doctype for automatic API
* **Quality code**: Correct project quality code errors with pre-commit automatic workflow (see below)
* **Factorization**: Factorize the code as much as possible. The command line tools shall only include the main workflow and rely on the cars python modules.
* **Be careful with user interface upgrade:** If major modifications of the user interface or of the tool's behaviour are done, update the user documentation (and the notebooks if necessary).
* **Logging and no print**: The usage of the `print()` function is forbidden: use the `logging` python standard module instead.
* **Limit classes**: If possible, limit the use of classes as much as possible and opt for a functional approach. The classes are reserved for data modelling if it is impossible to do so using `xarray` and for the good level of modularity.
* **Limit new dependencies**: Do not add new dependencies unless it is absolutely necessary, and only if it has a **permissive license**.

Pre-commit validation
=====================

A pre-commit validation is installed with code quality tools (see below).
It is installed automatically by `make install-dev` command.

Here is the way to install it manually:

.. code-block:: console

  $ pre-commit install

This installs the pre-commit hook in `.git/hooks/pre-commit`  from `.pre-commit-config.yaml <https://raw.githubusercontent.com/CNES/cars/master/.pre-commit-config.yaml>`_ file configuration.

It is possible to test pre-commit before commiting:

.. code-block:: console

  $ pre-commit run --all-files                # Run all hooks on all files
  $ pre-commit run --files shareloc/__init__.py   # Run all hooks on one file
  $ pre-commit run pylint                     # Run only pylint hook


Code quality
=============
Shareloc uses `Black`_,  and `Pylint`_ quality code checking.

Use the following command in CARS `virtualenv`_ to check the code with these tools:

.. code-block:: console

    $ make lint

Use the following command to format the code with black:

.. code-block:: console

    $ make format

Black
-----
`Black`_ is a quick and deterministic code formatter to help focus on the content.

Shareloc ``black`` configuration is done in `pyproject.toml`_

If necessary, Black doesn’t reformat blocks that start with "# fmt: off" and end with # fmt: on, or lines that ends with "# fmt: skip". "# fmt: on/off" have to be on the same level of indentation.

`Black`_ manual usage examples:

.. code-block:: console

    $ cd SHARELOC_HOME
    $ black --check shareloc tests  # Check code with black with no modifications
    $ black --diff shareloc tests   # Show black diff modifications
    $ black cars tests          # Apply modifications

Pylint
------
`Pylint`_ is a global linting tool which helps to have many information on source code.

CARS ``pylint`` configuration is done in dedicated `.pylintrc <http://https://raw.githubusercontent.com/CNES/cars/master/.pylintrc>`_ file.

`Pylint`_ messages can be avoided (in particular cases !) adding "# pylint: disable=error-message-name" in the file or line.
Look at examples in source code.

Pylint manual usage examples:

.. code-block:: console

  $ cd SHARELOC_HOME
  $ pylint tests shareloc       # Run all pylint tests
  $ pylint --list-msgs          # Get pylint detailed errors informations


Tests
======

Shareloc includes a set of tests executed with `pytest <https://docs.pytest.org/>`_ tool.

To launch tests:

.. code-block:: console

    make test

Advanced testing
----------------

To execute the tests manually, use ``pytest`` at the Shareloc projects's root (after initializing the environment):

.. code-block:: console

    $ python -m pytest

To run only the unit tests:

.. code-block:: console

    $ cd shareloc/
    $ pytest -m unit_tests



It is possible to obtain the code coverage level of the tests by installing the ``pytest-cov`` module and use the ``--cov`` option.

.. code-block:: console

    $ cd shareloc/
    $ python -m pytest --cov=cars

It is also possible to execute only a specific part of the test, either by indicating the test file to run:

.. code-block:: console

    $ cd shareloc/
    $ python -m pytest tests/test_tiling.py

Or by using the ``-k`` option which will execute the tests which names contain the option's value:

.. code-block:: console

    $ cd shareloc/
    $ python -m pytest -k triangulation

By default, ``pytest`` does not display the traces generated by the tests but only the tests' status (passed or failed). To get all traces, the following options have to be added to the command line (which can be combined with the previous options):

.. code-block:: console

    $ cd shareloc/
    $ python -m pytest -s -o log_cli=true -o log_cli_level=INFO


.. _`the GitHub repository`: https://github.com/CNES/shareloc
.. _`Shareloc Contribution guide`: https://github.com/CNES/Shareloc/blob/master/CONTRIBUTING.md
.. _`virtualenv`: https://virtualenv.pypa.io/
.. _`Black`: https://black.readthedocs.io/
.. _`Pylint`: http://pylint.pycqa.org/
.. _`pyproject.toml`: https://raw.githubusercontent.com/CNES/shareloc/master/pyproject.toml