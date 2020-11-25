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
import m2r
from recommonmark.parser import CommonMarkParser

sys.path.insert(0, os.path.abspath('../..'))
sys.setrecursionlimit(1500)

# -- Project information -----------------------------------------------------

project = 'Tree Box Model'
copyright = '2020, Olli-Pekka Tikkasalo'
author = 'Olli-Pekka Tikkasalo'

# The full version, including alpha/beta/rc tags
release = '0.3'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'recommonmark', 'sphinx.ext.todo',
              'sphinx.ext.mathjax']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# autodoc options
autodoc_member_order = 'alphabetical'
autodoc_parameter_order = 'alphabetical'

# napoleon options
napoleon_google_docstring = True
napoleon_use_rtype = False
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'nature'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# latex elements
latex_elements = {'preamble': r'''
\usepackage[T1]{fontenc}
\usepackage{textcomp}
\usepackage{lmodern}
\usepackage[utf8]{inputenc}
\DeclareUnicodeCharacter{03C1}{$\rho$}
\DeclareUnicodeCharacter{03C3}{$\sigma$}
\DeclareUnicodeCharacter{03B7}{$\eta$}
\DeclareUnicodeCharacter{03F5}{$\varepsilon$}
'''}


# Source suffices

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

source_parsers = {
    '.md': CommonMarkParser,
}


def docstring(app, what, name, obj, options, lines):
    md = '\n'.join(lines)
    rst = m2r.convert(md)
    lines.clear()
    lines += rst.splitlines()


def setup(app):
    app.connect('autodoc-process-docstring', docstring)
