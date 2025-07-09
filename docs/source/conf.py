from __future__ import annotations

import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

from sphinx.application import Sphinx

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))

import protdata


# -- General configuration ------------------------------------------------

project = "protdata"
author = "Max Frank"
copyright = f"{datetime.now():%Y}, Max Frank"
release = version = metadata.version("protdata")

templates_path = ["_templates"]
html_static_path = ["_static"]
source_suffix = {".rst": "restructuredtext", ".md": "myst-nb"}
master_doc = "index"
default_role = "literal"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
pygments_style = "sphinx"

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",
    "sphinx_issues",
    "sphinx_design",
    "sphinxext.opengraph",
    "scanpydoc",
    "IPython.sphinxext.ipython_console_highlighting",
]

myst_enable_extensions = ["colon_fence", "dollarmath"]
myst_heading_anchors = 3
nb_execution_mode = "off"

autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = True
napoleon_use_param = True
typehints_defaults = "braces"
issues_github_path = "czbiohub-sf/protdata"

nitpicky = True
nitpick_ignore = []

intersphinx_mapping = dict(
    python=("https://docs.python.org/3", None),
    numpy=("https://numpy.org/doc/stable", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable", None),
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    scipy=("https://docs.scipy.org/doc/scipy", None),
    h5py=("https://docs.h5py.org/en/latest", None),
)

# qualname_overrides = {
#     "h5py._hl.group.Group": "h5py.Group",
#     "h5py._hl.files.File": "h5py.File",
#     "h5py._hl.dataset.Dataset": "h5py.Dataset",
#     "anndata._core.anndata.AnnData": "anndata.AnnData",
# }


# -- Social cards ---------------------------------------------------------

# ogp_site_url = "https://protdata.readthedocs.io/"


# -- Options for HTML output ----------------------------------------------

html_theme = "scanpydoc"
html_theme_options = dict(
    use_repository_button=True,
    repository_url="https://github.com/czbiohub-sf/protdata",
    repository_branch="main",
    navigation_with_keys=False,
)
html_show_sphinx = False
html_logo = "_static/img/protdata_schema.svg"


def setup(app: Sphinx):
    app.add_css_file("css/custom.css")
