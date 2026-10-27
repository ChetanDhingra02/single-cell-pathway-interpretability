"""Compatibility shims so the site runs on older Streamlit releases.

The site is built and tested on Streamlit 1.66. Some keyword arguments it uses
(accessible `alt=` names, `required=` on segmented controls, and
`unsafe_allow_javascript=` on st.html) only exist in recent releases. These
wrappers pass them only when the installed version supports them. On older
versions the page still renders; the progressive-enhancement script (scroll
effects, interactive hero figure) simply does not run.
"""

from __future__ import annotations

import inspect

import streamlit as st


def _accepts(fn, name: str) -> bool:
    try:
        return name in inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False


HTML_JS = _accepts(st.html, "unsafe_allow_javascript")
PLOTLY_ALT = _accepts(st.plotly_chart, "alt")
DF_ALT = _accepts(st.dataframe, "alt")
SEG_REQUIRED = _accepts(st.segmented_control, "required")
TESTED = "1.66"


def html(body: str, js: bool = False) -> None:
    if js and HTML_JS:
        st.html(body, unsafe_allow_javascript=True)
    else:
        st.html(body)


def plotly(fig, alt: str | None = None, **kw) -> None:
    if alt and PLOTLY_ALT:
        kw["alt"] = alt
    st.plotly_chart(fig, **kw)


def dataframe(data, alt: str | None = None, **kw) -> None:
    if alt and DF_ALT:
        kw["alt"] = alt
    st.dataframe(data, **kw)


def seg(label: str, options, default, key: str, **kw):
    """Single-choice segmented control that always returns a value."""
    if SEG_REQUIRED:
        kw["required"] = True
    value = st.segmented_control(label, options, default=default, key=key, **kw)
    return default if value is None else value


def version_ok() -> bool:
    try:
        major, minor = (int(x) for x in st.__version__.split(".")[:2])
        return (major, minor) >= tuple(int(x) for x in TESTED.split("."))
    except ValueError:
        return True


def iframe(html_doc: str, height: int = 900) -> None:
    """Render a full HTML document in an iframe (st.iframe when available)."""
    if hasattr(st, "iframe"):
        st.iframe(html_doc, height=height)
    else:  # older Streamlit
        import streamlit.components.v1 as components

        components.html(html_doc, height=height, scrolling=False)
