Shareloc Documentation generation
=============================

Shareloc documentation is based on [Sphinx](https://www.sphinx-doc.org/) and is in `docs/source` directory.

Use the following command line at Shareloc source code root directly

```
make doc
```

for automated Shareloc installation and documentation generation in `docs/build`

Otherwise follow the steps below:


Shareloc installation with Sphinx dependencies
------------------------------------------

First, create a virtualenv and install Shareloc following [Shareloc Installation](./docs/source/install.rst)

You can use the following command line at CARS root directory:

```
make install-doc
```

This installs Shareloc in a virtualenv with sphinx documentation dependencies using : `pip install .[doc]`  

The autodoc needs indeed Shareloc installed for modules API.


Shareloc Sphinx documentation
-------------------------

Go to `docs` directory from Shareloc source root directory.

```
cd docs
```

First generate Autodoc generation in docs directory:
````
sphinx-apidoc -o source/apidoc/ ../shareloc
````

For HTML documentation generation:
```
make html
```

For PDF generation :
```
make latexpdf
```

To clean generated documentation :
```
make clean
```

Plantuml diagrams
-----------------

To generate/update CARS [plantUML](https://plantuml.com/) documentation diagrams:
* Install plantuml on your system.
* Go to diagram directory
* Use following command line to generate in SVG:
```
plantuml -tsvg diagram.pu
```

Diagrams source and associated SVG are in `docs/source/diagrams/` directory
