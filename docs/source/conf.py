# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))  # adjust path if needed


project = 'ideal-genom'
copyright = '2025, Luis Giraldo Gonzalez Ricardo'
author = 'Luis Giraldo Gonzalez Ricardo'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
"sphinx.ext.autodoc",
"sphinx.ext.napoleon", # for NumPy/Google-style docstrings\
"sphinx_autodoc_typehints",
"myst_parser", # if you use Markdown\
]

autodoc_typehints = "description"


templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
