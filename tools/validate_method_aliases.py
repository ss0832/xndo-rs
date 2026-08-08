#!/usr/bin/env python3
"""Static release guard for historically distinct INDO vs INDO/S aliases."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
text = (root / "src" / "method.rs").read_text(encoding="utf-8")

# The bare MolDS label must resolve to the separate ground-state method.
assert '"indo" => Some(Self::Indo)' in text
# Explicit INDO/S spellings may resolve to ZINDO/S, but bare INDO may not.
z_accept = re.search(r'Self::ZindoS\s*=>\s*&\[(.*?)\]', text, re.S)
assert z_accept, "ZINDO/S accepted_strings entry missing"
assert '"indo"' not in z_accept.group(1), "bare indo leaked into ZINDO/S accepted_strings"
assert '"indo/s"' in z_accept.group(1)
assert '"indos"' in z_accept.group(1)
indo_accept = re.search(r'Self::Indo\s*=>\s*&\[(.*?)\]', text, re.S)
assert indo_accept and '"indo"' in indo_accept.group(1)

for rel in ["README.md", "docs/python-api.md", "docs/rust-api.md"]:
    doc = (root / rel).read_text(encoding="utf-8")
    # Every executable ZINDO/S alias row must avoid a comma-delimited bare indo.
    for line in doc.splitlines():
        if "ZINDO/S" in line and "Accepted" not in line and "`zindo/s`" in line:
            assert "`indo`," not in line and ", `indo`" not in line, f"bare indo advertised as ZINDO/S in {rel}: {line}"

print("method alias audit: OK (indo != ZINDO/S; indo/s remains explicit INDO/S alias)")
