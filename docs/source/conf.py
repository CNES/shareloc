# pylint: skip-file
# flake8: noqa
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

# Extend Recursion limit for RecursionError in big files (bug astroid)
sys.setrecursionlimit(8 * sys.getrecursionlimit())

# -- Project information -----------------------------------------------------

project = "Shareloc"
copyright = "2023, CNES"
author = "Shareloc Team"

# The full version, including alpha/beta/rc tags
from pkg_resources import get_distribution

try:
    version = get_distribution("Shareloc").version
    release = version
except Exception as error:
    print("WARNING: cannot find shareloc version")
    version = "Unknown"
    release = version

# The master toctree document.
master_doc = "index"

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.ifconfig",  # add if config possibility in rst files
    "sphinx.ext.intersphinx",  # other projects automatic links to doc
    "sphinx.ext.imgmath",
    "sphinx.ext.autodoc",  # apidoc automatic generation
    "sphinx.ext.viewcode",  # viewcode in automatic apidoc
    "sphinx.ext.autosectionlabel",  # label automatic in sections.
    "autoapi.extension",
    "sphinx_tabs.tabs",
]

# imgmath configuration
imgmath_embed = True

# Autoapi configuration
autoapi_dirs = ["../../shareloc"]
autoapi_root = "api_reference"
autoapi_keep_files = True
autoapi_options = [
    "members",
    "undoc-members",
    "private-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
]

# add autohint
autodoc_typehints = "description"

# remove warning autosectionlabel
suppress_warnings = ["autosectionlabel.*"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Title
html_title = "Shareloc Documentation"
html_short_title = "Shareloc Documentation"

# Logo
html_logo = "images/shareloc_picto.png"

# Favicon
html_favicon = "images/favicon_shareloc.ico"

# Theme options
html_theme_options = {
    "logo_only": True,
    "navigation_depth": 3,
}

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
# html_css_files = ["css/my_custom.css"]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "SharelocDoc"


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    "papersize": "letterpaper",
    # The font size ('10pt', '11pt' or '12pt').
    "pointsize": "10pt",
    # Additional stuff for the LaTeX preamble.
    "preamble": "",
    # Latex figure (float) alignment
    "figure_align": "htbp",
}
numfig = True

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "SharelocDoc.tex", "Shareloc documentation", "Rich Yap", "manual"),
]

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "sharelocdocs", "Shareloc Documentation", [author], 1)]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "SharelocDoc",
        "Shareloc Documentation",
        author,
        "SharelocDoc",
        "One line description of project.",
        "Miscellaneous",
    ),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ["search.html"]
