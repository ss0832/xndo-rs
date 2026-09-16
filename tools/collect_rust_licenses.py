#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Collect the licence texts of every Rust crate that ends up in a build.

Writes `THIRD_PARTY_LICENSES_RUST.md`. Run with `--check` in CI to fail when it
is stale.

**Why this is not `cargo-about`.** Generators of that kind synthesise licence
text from an SPDX store keyed on the crate's declared `license` field, rather
than reading what the crate actually ships. `faer` is the case that settles it
here: its manifest says `MIT`, but the crate bundles

    COPYING.EIGEN.MPL2
    COPYING.LAPACK.BSD
    COPYING.SUITE_SPARSE.AMD.BSD
    COPYING.SUITE_SPARSE.COLAMD.BSD

because parts of it derive from Eigen, LAPACK and SuiteSparse. A tool that
trusts the SPDX field writes "MIT" and **drops the MPL-2.0 obligation entirely**.
`faer` is a direct dependency and is in every binary this project produces, so
that is not a hypothetical.

This script instead walks the unpacked crate sources in the Cargo registry and
copies `LICENSE*`, `LICENCE*`, `COPYING*`, `NOTICE*` and `AUTHORS*` verbatim. It
needs no network and no tool beyond the `cargo metadata` that ships with cargo.

Linkage is classified, because it changes what has to be conveyed: a crate
linked into the library or the binary is distributed with it, while a
`proc-macro` or build-dependency runs at build time and is not. Both are listed
-- the distinction is recorded, not used to hide anything -- because `pyo3`
pulls a proc-macro subtree that would otherwise look like shipped code.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "THIRD_PARTY_LICENSES_RUST.md")

LICENSE_FILE = re.compile(r"^(LICEN[CS]E|COPYING|NOTICE|AUTHORS)", re.I)
# A licence file that is plainly a duplicate of another in the same crate adds
# nothing; the dual-licence pair is kept because both halves are offered.
MAX_BYTES = 200_000


def cargo_metadata():
    out = subprocess.run(
        ["cargo", "metadata", "--locked", "--format-version", "1"],
        capture_output=True, cwd=ROOT, check=True,
    ).stdout
    return json.loads(out.decode("utf-8-sig"))


def classify_linkage(meta):
    """Split packages into what is linked into the product and what only builds it.

    `proc-macro` crates and anything reachable only through a build-dependency
    run on the build machine. Everything else reachable from the root package's
    normal dependencies is linked in and therefore conveyed.
    """
    packages = {p["id"]: p for p in meta["packages"]}
    nodes = {n["id"]: n for n in meta["resolve"]["nodes"]}
    roots = [m for m in meta["workspace_members"]]

    linked, buildtime = set(), set()

    def walk(pid, kinds, seen):
        if pid in seen:
            return
        seen.add(pid)
        pkg = packages[pid]
        is_macro = any(t["kind"] == ["proc-macro"] for t in pkg.get("targets", []))
        if "build" in kinds or is_macro:
            buildtime.add(pid)
            next_kinds = kinds | {"build"}
        else:
            linked.add(pid)
            next_kinds = kinds
        for dep in nodes[pid]["deps"]:
            for dk in dep["dep_kinds"]:
                kind = dk.get("kind")
                if kind == "dev":
                    continue
                walk(dep["pkg"], next_kinds | ({"build"} if kind == "build" else set()), seen)

    for root in roots:
        walk(root, set(), set())
    # A crate reachable both ways is conveyed, so linkage wins.
    buildtime -= linked
    return packages, linked, buildtime


def license_files(pkg):
    """Every licence-ish file the crate actually ships, read verbatim."""
    manifest = pkg.get("manifest_path")
    if not manifest:
        return []
    src = os.path.dirname(manifest)
    found = []
    try:
        entries = sorted(os.listdir(src))
    except OSError:
        return []
    for entry in entries:
        if not LICENSE_FILE.match(entry):
            continue
        path = os.path.join(src, entry)
        if not os.path.isfile(path) or os.path.getsize(path) > MAX_BYTES:
            continue
        with io.open(path, encoding="utf-8", errors="replace") as fh:
            found.append((entry, fh.read()))
    return found


def render(meta, packages, linked, buildtime, lock_sha):
    """Render the notices, with each distinct licence text reproduced once.

    Forty-nine crates ship a byte-identical copy of Apache-2.0. Reproducing it
    forty-nine times made this file 864 KB and unreadable, and conveyed nothing
    that one copy plus the list of crates does not. Texts are therefore grouped
    by content hash: every crate is still named, every distinct text is still
    verbatim, and a crate whose copy differs by so much as a copyright line gets
    its own entry.
    """
    root_names = {packages[m]["name"] for m in meta["workspace_members"]}

    # content hash -> (text, [(crate, version, filename)])
    texts = {}
    rows = []
    for group, linkage in ((linked, "linked"), (buildtime, "build-time")):
        # Sorted: iterating a set walks it in an order that depends on the
        # per-process string hash seed, which made this file differ between runs
        # and `--check` fail against output that was not wrong.
        for pid in sorted(group, key=lambda i: (packages[i]["name"], packages[i]["version"])):
            pkg = packages[pid]
            if pkg["name"] in root_names:
                continue
            files = license_files(pkg)
            digests = []
            for filename, text in files:
                digest = hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]
                entry = texts.setdefault(digest, (text, []))
                entry[1].append((pkg["name"], pkg["version"], filename))
                digests.append((filename, digest))
            rows.append((pkg["name"], pkg["version"], linkage,
                         pkg.get("license") or "(none declared)",
                         pkg.get("repository") or "", digests))
    rows.sort()

    out = []
    out.append("# Third-party Rust crate licences")
    out.append("")
    out.append("<!-- cargo-lock-sha256: {} -->".format(lock_sha))
    out.append("")
    out.append(
        "Generated by `tools/collect_rust_licenses.py` from the unpacked crate\n"
        "sources in the Cargo registry. The texts are **verbatim copies of the\n"
        "files each crate ships**, not text synthesised from its declared SPDX\n"
        "expression -- see the script's header for why that distinction matters\n"
        "for `faer`."
    )
    out.append("")
    out.append(
        "Each distinct licence text appears once, in the second half, with every\n"
        "crate that ships it named. A crate whose copy differs by so much as a\n"
        "copyright line gets its own entry."
    )
    out.append("")
    out.append(
        "Regenerate with `python tools/collect_rust_licenses.py`; CI runs it with\n"
        "`--check` so a dependency change not reflected here fails."
    )
    out.append("")
    out.append("## Crates")
    out.append("")
    out.append(
        "`linked` crates are conveyed with every artifact this project produces.\n"
        "`build-time` crates are proc-macros and build dependencies: they run on\n"
        "the build machine and are not part of the distributed artifacts. Both are\n"
        "listed; the distinction is recorded, not used to omit anything."
    )
    out.append("")
    out.append("| crate | version | linkage | declared | licence files |")
    out.append("| --- | --- | --- | --- | --- |")
    for name, version, linkage, declared, repo, digests in rows:
        if digests:
            files = ", ".join(
                "{} ([{}](#text-{}))".format(f, d, d) for f, d in digests)
        else:
            files = "**none shipped**; the declared expression is the only statement"
        out.append("| `{}` | {} | {} | `{}` | {} |".format(
            name, version, linkage, declared, files))
    out.append("")
    out.append("## Licence texts")
    out.append("")
    for digest, (text, users) in sorted(texts.items(), key=lambda kv: -len(kv[1][1])):
        crates = sorted({"{} {}".format(n, v) for n, v, _ in users})
        filenames = sorted({f for _, _, f in users})
        out.append("### text-{}".format(digest))
        out.append("")
        out.append("Shipped as {} by {} crate(s): {}".format(
            ", ".join("`{}`".format(f) for f in filenames),
            len(crates), ", ".join("`{}`".format(c) for c in crates)))
        out.append("")
        out.append("```")
        out.append(text.rstrip("\n"))
        out.append("```")
        out.append("")
    return "\n".join(out).rstrip("\n") + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="regenerate and fail if the committed file differs")
    args = ap.parse_args()

    with open(os.path.join(ROOT, "Cargo.lock"), "rb") as fh:
        lock_sha = hashlib.sha256(fh.read()).hexdigest()

    meta = cargo_metadata()
    packages, linked, buildtime = classify_linkage(meta)
    text = render(meta, packages, linked, buildtime, lock_sha)

    if args.check:
        existing = ""
        if os.path.exists(OUT):
            with io.open(OUT, encoding="utf-8") as fh:
                existing = fh.read()
        if existing != text:
            print(
                "THIRD_PARTY_LICENSES_RUST.md is stale; run "
                "`python tools/collect_rust_licenses.py`",
                file=sys.stderr,
            )
            return 1
        print("THIRD_PARTY_LICENSES_RUST.md is up to date "
              "({} linked, {} build-time)".format(len(linked) - 1, len(buildtime)))
        return 0

    with io.open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("wrote {} ({} linked, {} build-time)".format(
        os.path.relpath(OUT, ROOT), len(linked) - 1, len(buildtime)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
