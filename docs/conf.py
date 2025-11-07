# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath("../jayrock"))


# -- Project information -----------------------------------------------------

project = "jayrock"
copyright = "2025, Max Mahlke"
author = "Max Mahlke"

# The full version, including alpha/beta/rc tags
release = "0.1"
html_title = "jayrock"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_prompt",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
    "sphinx.ext.graphviz",
    "hoverxref.extension",
    "sphinx_design",
    "sphinx_copybutton",
    "myst_parser",
    "sphinx-togglebutton",
]

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]
myst_url_schemes = ("http", "https", "mailto")

# Print out todos in documentation?
todo_include_todos = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


graphviz_output_format = "svg"


# ------
# hoverxref_role_types = {
#     "hoverxref": "modal",
#     "ref": "modal",  # for hoverxref_auto_ref config
#     "confval": "tooltip",  # for custom object
#     "mod": "tooltip",  # for Python Sphinx Domain
#     "class": "tooltip",  # for Python Sphinx Domain
#     "term": "tooltip",  # for Python Sphinx Domain
# }

hoverxref_roles = ["numref", "confval", "setting", "term"]

hoverxref_project = "rocks"
hoverxref_version = "latest"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
# html_logo = "_static/logo_jayrock.svg"
html_title = "jayrock Documentation"

pygments_style = "nord"
pygments_dark_style = "nord-darker"

html_show_copyright = False
html_show_sphinx = False

html_theme_options = {
    "logo": {
        "image_light": "_static/logo_jayrock.svg",
        "image_dark": "_static/logo_jayrock_dark.svg",
    },
    # "light_logo": "logo_jayrock.svg",
    # "dark_logo": "logo_jayrock_dark.svg",
    "repository_url": "https://github.com/maxmahlke/jayrock",
    "use_repository_button": True,
    "use_sidenotes": False,
    "show_toc_level": 3,
}

html_last_upadted_fmt = "%b %d, %Y"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_style = "css/custom.css"

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    "css/custom.css",
]
