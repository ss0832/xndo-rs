#!/usr/bin/env python3
"""Static release guard for CNDO/2 and MolDS-labelled INDO wiring.

This cannot replace `cargo test`; it prevents API/dispatch/docs drift on hosts
where a Rust toolchain is unavailable.
"""
from pathlib import Path
import csv
import re

ROOT = Path(__file__).resolve().parents[1]

def text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")

method = text("src/method.rs")
engine = text("src/cndo_indo.rs")
xndo = text("src/xndo.rs")
python = text("src/python.rs")
cli = text("src/bin/xndo_rs.rs")
methods_doc = text("docs/methods.md")
python_doc = text("docs/python-api.md")

# Native method registry and aliases.
assert "Self::Cndo2 | Self::Indo | Self::Mndo" in method
assert 'Self::Cndo2 => &["cndo2", "cndo/2", "cndo"]' in method
assert 'Self::Indo => &["indo"]' in method
assert '"indo" => Some(Self::Indo)' in method
assert '"indos" | "indo/s" | "zindo" | "zindos" | "zindo/s" => Some(Self::ZindoS)' in method
ground_api = re.search(
    r'Self::Cndo2 \| Self::Indo \| Self::Mindo3 => &\[(.*?)\]', method, re.S
)
assert ground_api, "CNDO/2 and INDO API registry entry missing"
for api in ["single_point", "gradient", "forces", "hessian", "frequencies"]:
    assert f'"{api}"' in ground_api.group(1)

# Real engine is public and independently dispatched at every exposed layer.
assert "pub fn run_cndo_indo(" in engine
assert "Method::Cndo2 => build_fock_cndo_rhf" in engine
assert "Method::Indo => build_fock_indo_rhf" in engine
assert "Method::Cndo2 => build_fock_cndo_uhf" in engine
assert "Method::Indo => build_fock_indo_uhf" in engine
assert "Method::Cndo2 | Method::Indo =>" in xndo
assert "CalculationResult::CndoIndo" in xndo
assert "Method::Cndo2 | Method::Indo =>" in python
assert "run_cndo_indo(&mol, meth" in python
assert "Method::Cndo2 | Method::Indo =>" in cli
assert "run_cndo_indo(&molecule, cli.method" in cli

# Runtime donor table scope is exactly H/Li/C/N/O/S; INDO itself gates sulfur.
raw = text("src/data/molds_cndo2_indo_parameters.csv")
rows = list(csv.DictReader([ln for ln in raw.splitlines() if ln and not ln.startswith("#")]))
assert [int(row["z"]) for row in rows] == [1, 3, 6, 7, 8, 16]
assert "if matches!(z, 1 | 3 | 6 | 7 | 8)" in engine
assert "MolDS-compatible INDO in xndo-rs 0.2.4 is enabled for H/Li/C/N/O" in engine

# An independent published H2 numeric oracle is frozen as an executable test.
assert "h2_matches_published_cndo2_indo_oracle" in engine
assert "1.474625" in engine

# User-visible documentation must state the method separation explicitly.
assert "| CNDO/2 | native | H, Li, C, N, O, S |" in methods_doc
assert "| INDO | native | H, Li, C, N, O |" in methods_doc
assert "not an\nalias for `INDO/S`" in methods_doc
assert "Bare `indo` selects the\nground-state INDO engine" in python_doc

print("CNDO/2 + INDO wiring audit: OK")
