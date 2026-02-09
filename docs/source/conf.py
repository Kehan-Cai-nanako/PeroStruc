# docs/source/conf.py
from __future__ import annotations

import os
import sys
from datetime import datetime

### Make your package importable for autodoc
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..", "src")))

### Project metadata shown in the built docs
project = "Disordered Perovskite Monte Carlo Optimizer"
author = "Alyssa (Xinyu) Xu, Kehan Cai, Pinchen Xie"
copyright = f"{datetime.now().year}, {author}"

### Pull version from your package. If you have __version__ in src/my_package/__init__.py, this will work.
try:
    import perostruc  # noqa: WPS433
    release = getattr(perostruc, "__version__", "0.1.0")
except Exception:
    release = "0.1.0"

### Extensions: add features to Sphinx
extensions = [
    "sphinx.ext.autodoc",       # generate docs from docstrings
    "sphinx.ext.napoleon",      # understand Google/NumPy docstring styles
    "sphinx.ext.autosummary",   # create summary tables + stub pages
    "sphinx.ext.viewcode",      # add links to highlighted source code
    "sphinx.ext.mathjax",       # render math formulas in HTML
    "myst_parser",              # allow Markdown (.md) pages if you want
]

### Templates / excludes
templates_path = ["_templates"]
exclude_patterns = []

### autodoc/autosummary behavior
autosummary_generate = True         # auto-create files under api "generated/"
autodoc_typehints = "description"   # show type hints in text, not in signatures
napoleon_google_docstring = True
napoleon_numpy_docstring = True

### HTML output options
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
