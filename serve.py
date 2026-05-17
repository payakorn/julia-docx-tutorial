#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["livereload"]
# ///
"""
Local dev server for the lecture site.

Run from the project root:
    ./serve.py            # http://localhost:8000/
    ./serve.py 9000       # custom port

Serves site/ as the document root, so the page lives at the site root
(/) rather than under /site/. Auto-reloads the browser when files
inside site/ change — including site/output/figures/, where the figure
generator now writes directly.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from livereload import Server

ROOT = Path(__file__).resolve().parent
SITE = ROOT / "site"
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8000


def write_version_json() -> None:
    """Write site/version.json with git SHA + server start time.

    Reads .git directly so we don't need the git binary in the image.
    """
    git_dir = ROOT / ".git"
    sha = "unknown"
    try:
        head = (git_dir / "HEAD").read_text().strip()
        if head.startswith("ref: "):
            ref = head[5:]
            sha = (git_dir / ref).read_text().strip()
        else:
            sha = head
    except OSError:
        pass

    payload = {
        "commit": sha,
        "commit_short": sha[:7] if sha != "unknown" else sha,
        "started_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (SITE / "version.json").write_text(json.dumps(payload) + "\n")


write_version_json()

server = Server()
server.watch(str(SITE))
print(f"→  http://localhost:{PORT}/   (Ctrl-C to stop)")
server.serve(root=str(SITE), port=PORT, host="0.0.0.0")
