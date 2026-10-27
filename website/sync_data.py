"""Refresh website/data_snapshot/ from the repository's frozen result files and
rebuild the shipped data bundle (data_snapshot/payload.json).

The site prefers the repository's own results/ when it runs inside the repo,
and falls back to this snapshot when the website folder is used on its own
(e.g. website-only.zip or a separate deployment). Run from anywhere:

    python website/sync_data.py
"""

import shutil
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent
REPO = APP.parent
DST = APP / "data_snapshot"
sys.path.insert(0, str(APP))

from ui.paths import ALL_FILES  # noqa: E402

missing = [f for f in ALL_FILES if not (REPO / f).exists()]
if missing:
    raise SystemExit(f"Not inside the repository (missing {missing[0]}); nothing copied.")
for f in ALL_FILES:
    (DST / f).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO / f, DST / f)
print(f"Copied {len(ALL_FILES)} files to {DST}")

from ui.paths import fingerprint  # noqa: E402
from ui.site import SHIPPED_BUNDLE, build_bundle, write_bundle  # noqa: E402

write_bundle(SHIPPED_BUNDLE, fingerprint(), build_bundle())
print(f"Rebuilt {SHIPPED_BUNDLE}")
