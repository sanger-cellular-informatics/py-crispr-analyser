# Copyright (C) 2025 Genome Research Ltd.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
from sphinx_pyproject import SphinxConfig

sys.path.insert(0, os.path.abspath('../../py_crispr_analyser/'))

config = SphinxConfig("../../pyproject.toml", globalns=globals())

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser', 'sphinx.ext.autodoc']
myst_enable_extensions = ['substitution']
templates_path = ['_templates']
exclude_patterns = []
autodoc_typehints = 'description'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
