# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'hy-vent'
copyright = '2026, Jonathan Mette'
author = 'Jonathan Mette'
release = '1.1.0'

# -- Point to python code ---------------------------------------------------

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '..', 'src', 'hy-vent').resolve()))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions =  [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',   
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#preinstalled fallback theme
#html_theme = "alabaster"
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
