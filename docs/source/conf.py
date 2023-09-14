# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "subproptools"
copyright = "2023, Kevin Lefrancois-Gagnon, Robert Mawhinney"
author = "Kevin Lefrancois-Gagnon, Robert Mawhinney"
release = "0.2.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ["_templates"]
exclude_patterns = []


import os

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
import pathlib
import sys

sys.path.insert(0, os.path.abspath("../.."))
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
html_theme = "alabaster"
html_static_path = ["_static"]
extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
]
