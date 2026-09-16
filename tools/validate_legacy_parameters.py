#!/usr/bin/env python3
"""Validate the provenance-audited parameter tables retained in xndo-rs.

The version is read from Cargo.toml rather than written here. A hardcoded
one goes stale at the first release and then quietly contradicts the tree.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "data"


def crate_version() -> str:
    """The version in Cargo.toml, which is the single source of truth."""
    for line in (ROOT / "Cargo.toml").read_text(encoding="utf-8").splitlines():
        if line.startswith("version"):
            return line.split("=", 1)[1].strip().strip('"')
    raise SystemExit("Cargo.toml has no version")

EXPECTED_ROWS = {
    "molds_cndo2_indo_parameters.csv": 6,
    "molds_zindo_s_parameters.csv": 5,
    "mindo3_parameters.csv": 10,
    "mindo3_pair_parameters.csv": 40,
    "legacy_parameter_catalog.csv": 7,
}

MANIFEST_FILES = [
    "element_data.csv",
    "mndo_parameters.csv",
    "mndo_pair_parameters.csv",
    "mndod_parameters.csv",
    "mndod_pair_parameters.csv",
    "zindo_s_parameters.csv",
    "molds_cndo2_indo_parameters.csv",
    "molds_zindo_s_parameters.csv",
    "mindo3_parameters.csv",
    "mindo3_pair_parameters.csv",
    "legacy_parameter_catalog.csv",
]


def rows(path: Path) -> list[dict[str, str]]:
    text = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return list(csv.DictReader(text))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    for name, count in EXPECTED_ROWS.items():
        got = rows(DATA / name)
        if len(got) != count:
            raise SystemExit(f"{name}: expected {count} data rows, got {len(got)}")

    pair_rows = rows(DATA / "mindo3_pair_parameters.csv")
    pair = {
        (int(row["z1"]), int(row["z2"])): (row["beta_ab"], row["alpha_ab"])
        for row in pair_rows
    }
    expected_tail = {
        (15, 17): ("0.277322", "1.543720"),
        (16, 17): ("0.221764", "1.950318"),
        (17, 17): ("0.258969", "1.792125"),
    }
    for key, expected in expected_tail.items():
        if pair.get(key) != expected:
            raise SystemExit(f"MINDO/3 pair {key}: expected {expected}, got {pair.get(key)}")

    catalog = {row["dataset_id"]: row for row in rows(DATA / "legacy_parameter_catalog.csv")}
    if catalog["molds_cndo2_indo"]["upstream_license"] != "GPL-3.0-or-later":
        raise SystemExit("MolDS CNDO2/INDO license marker changed")
    if catalog["mopac7_mindo3"]["upstream_license"] != "Public-Domain":
        raise SystemExit("MOPAC7 MINDO/3 license marker changed")

    cndo_indo = {
        row["sym"]: row
        for row in rows(DATA / "molds_cndo2_indo_parameters.csv")
    }
    lithium = cndo_indo.get("Li")
    if lithium is None or lithium["bonding_parameter_ev"] != "-9.0" or lithium["z_eff_l"] != "1.3":
        raise SystemExit("MolDS lithium parameter sentinel changed")

    manifest_path = DATA / "legacy_parameter_manifest.sha256"
    expected_manifest = {
        name: sha256(DATA / name) for name in MANIFEST_FILES
    }
    manifest = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        digest, name = line.split(None, 1)
        manifest[name.strip()] = digest
    if manifest != expected_manifest:
        missing = sorted(set(expected_manifest) - set(manifest))
        extra = sorted(set(manifest) - set(expected_manifest))
        changed = sorted(k for k in set(manifest) & set(expected_manifest) if manifest[k] != expected_manifest[k])
        raise SystemExit(f"manifest mismatch: missing={missing} extra={extra} changed={changed}")

    print("legacy parameter audit: OK")
    for name in MANIFEST_FILES:
        print(f"  {name}: {expected_manifest[name]}")


if __name__ == "__main__":
    main()
