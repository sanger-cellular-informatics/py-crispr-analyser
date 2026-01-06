# Copyright (C) 2025 Genome Research Ltd.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
from sphinx_pyproject import SphinxConfig

# Put project root (so package can be imported as `py_crispr_analyser`) on sys.path
sys.path.insert(0, os.path.abspath('../../'))

config = SphinxConfig("../../pyproject.toml", globalns=globals())

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Enable autosummary (generates pages with signatures), napoleon for Google/NumPy docstrings,
# and viewcode to link to source, in addition to the existing autodoc.
extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

# Generate autosummary pages automatically
autosummary_generate = True

# If you have heavy optional dependencies that are not available when building docs (e.g. on RTD),
# mock them so modules can be imported for autodoc.
autodoc_mock_imports = [
    'numba',
    'numba.cuda',
    'numpy',
]

# Default autodoc options so that members and undocumented members are included
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

myst_enable_extensions = ['substitution']
templates_path = ['_templates']
exclude_patterns = []
autodoc_typehints = 'description'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
