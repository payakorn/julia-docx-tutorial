#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["livereload"]
# ///
"""
Local dev server for the lecture site.

Run from the project root:
    ./serve.py            # http://localhost:8000/site/
    ./serve.py 9000       # custom port

Auto-reloads the browser when site/ or output/figures/ changes.
"""
import sys
from pathlib import Path
from livereload import Server

ROOT = Path(__file__).resolve().parent
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8000

server = Server()
server.watch(str(ROOT / "site"))
server.watch(str(ROOT / "output" / "figures"))
print(f"→  http://localhost:{PORT}/site/   (Ctrl-C to stop)")
server.serve(root=str(ROOT), port=PORT, host="0.0.0.0")
