"""Assembles the single-page front end (web/) with the data bundle and shows it
full-screen inside Streamlit.

Startup is fast because the data bundle is stored on disk as JSON together with
a fingerprint of every input table and of the bundle-building code
(ui/paths.py). When the fingerprint matches, pandas is never imported. When a
result file or the bundle code changes, the bundle is rebuilt once from the
frozen tables (ui/payload.py) and saved again.
"""

from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from .compat import iframe
from .paths import APP, fingerprint

WEB = APP / "web"
LOCAL_BUNDLE = APP / ".cache" / "payload.json"         # written at runtime (git-ignored)
SHIPPED_BUNDLE = APP / "data_snapshot" / "payload.json"  # shipped with the site, refreshed by sync_data.py

HOST_CSS = """<style>
[data-testid="stHeader"], [data-testid="stToolbar"], [data-testid="stDecoration"], [data-testid="stStatusWidget"], footer { display:none !important; }
html, body, [data-testid="stApp"], [data-testid="stAppViewContainer"], [data-testid="stMain"] { overflow:hidden !important; background:#fbf8f1 !important; }
[data-testid="stMainBlockContainer"] { padding:0 !important; max-width:none !important; }
[data-testid="stVerticalBlock"] { gap:0 !important; }
iframe { position:fixed !important; inset:0 !important; width:100vw !important; height:100vh !important; height:100dvh !important;
         border:0 !important; z-index:999999; display:block; background:#fbf8f1; }
</style>"""


def _read_bundle(path: Path, fp: str) -> str | None:
    try:
        head, _, body = path.read_text(encoding="utf-8").partition("\n")
    except OSError:
        return None
    return body if head.strip() == fp and body else None


def build_bundle() -> str:
    from .payload import build_payload  # imports pandas only when a rebuild is needed

    return json.dumps(build_payload(), separators=(",", ":"), ensure_ascii=False)


def write_bundle(path: Path, fp: str, data: str) -> bool:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(fp + "\n" + data, encoding="utf-8")
        return True
    except OSError:  # read-only deployment: just keep it in memory
        return False


@st.cache_resource(show_spinner=False)
def bundle_json() -> str:
    fp = fingerprint()
    for path in (LOCAL_BUNDLE, SHIPPED_BUNDLE):
        data = _read_bundle(path, fp)
        if data is not None:
            return data
    data = build_bundle()
    write_bundle(LOCAL_BUNDLE, fp, data)
    return data


def _web_stamp() -> tuple:
    return tuple((p.name, p.stat().st_mtime_ns) for p in sorted(WEB.iterdir()) if p.is_file())


@st.cache_resource(show_spinner=False)
def page_html(stamp: tuple) -> str:
    html = (WEB / "index.html").read_text(encoding="utf-8")
    css = (WEB / "site.css").read_text(encoding="utf-8")
    js = (WEB / "site.js").read_text(encoding="utf-8")
    data = bundle_json().replace("</", "<\\/")
    return html.replace("/*__CSS__*/", css).replace("/*__DATA__*/null", data).replace("/*__JS__*/", js)


def build_html() -> str:
    return page_html(_web_stamp())


def render_site() -> None:
    st.html(HOST_CSS)
    iframe(build_html(), height=900)
