"""Research website: Auditing Faithfulness and Overlap in Pathway-Constrained Single-Cell Models.

Run from the repository root (so .streamlit/config.toml is picked up):

    pip install -r website/requirements.txt
    streamlit run website/streamlit_app.py

Python reads the frozen result tables (results/ in the repo, or the bundled
copy in website/data_snapshot/), packs them into one JSON bundle (cached on disk
and keyed by a fingerprint of the inputs) and serves the interactive front end
in website/web/. No raw data or model weights are loaded.
"""

import sys
from pathlib import Path

import streamlit as st

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

st.set_page_config(page_title="Pathway audit", page_icon=":material/genetics:", layout="wide",
                   initial_sidebar_state="collapsed")

from ui.paths import missing_files  # noqa: E402
from ui.site import render_site  # noqa: E402

missing = missing_files()
if missing:
    st.error(f"A result file the site needs is missing: `{missing[0]}`. Keep `results/` and `data/gene_sets/` next to "
             "`website/`, or the bundled `website/data_snapshot/` inside it.", icon=":material/error:")
    st.stop()

render_site()
