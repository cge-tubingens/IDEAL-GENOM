import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

project = 'ideal-genom'
copyright = '2025, Luis Giraldo Gonzalez Ricardo'
author = 'Luis Giraldo Gonzalez Ricardo'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    #'no-index': True,  # Avoid duplicate warnings
}

autodoc_class_signature = 'separated'

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = False

templates_path = ['_templates']
exclude_patterns = []

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
