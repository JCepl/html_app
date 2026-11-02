#!/usr/bin/env python3
"""Cache-busting: put a short content hash on app.css / app.js in index.html (?v=...).

Run after editing app.css or app.js and before committing:  python3 stamp_version.py
The browser then re-downloads a file only when its content actually changed.
"""
import hashlib, re, pathlib

root = pathlib.Path(__file__).parent
html = (root / "index.html").read_text()
for name in ("app.css", "app.js"):
    h = hashlib.sha1((root / name).read_bytes()).hexdigest()[:10]
    html = re.sub(r'(%s)\?v=[0-9a-f]+' % re.escape(name), r'\1?v=' + h, html)
(root / "index.html").write_text(html)
print("stamped")
