#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Scan the tree, or built artifacts, for personal information.

Absolute paths that contain a user name, a home directory, a machine name or an
account identifier must never be written into anything that is committed or
shipped. The one people miss is **generated binary artifacts**: a compiler or a
build tool that stamps its own command line into an output file embeds whatever
path it was given, so the artifact has to be checked and not just the source.

    python tools/privacy_scan.py --tree .
    python tools/privacy_scan.py --artifacts dist/*.whl target/package/*.crate ...

Archives are read **without extracting** (`.zip`, `.whl`, `.crate`, `.tar.gz`),
and binaries are searched in both UTF-8 and UTF-16LE, because a Windows
toolchain writes paths in the latter and a UTF-8-only scan reports a clean bill
of health on a file that is full of them.

**The identity patterns are resolved at run time** from `getpass.getuser()`,
`socket.gethostname()` and the environment, so this file contains no personal
information of its own and works on any machine. Reports go to stdout only: a
scanner that writes a report containing full paths back into the repository
would be defeating itself.

## What is allowed, and where

Allowlisting is per **document**, not per pattern. Upstream author names and
e-mail addresses inside licence and attribution files are *required* to be
retained -- removing them breaks Apache-2.0 4(c) and the GPL -- so the `email`
class is permitted inside those files and nowhere else. Path classes are
permitted nowhere at all, and the `local-identity` class is permitted nowhere,
because that is the one that is actually a privacy risk rather than a tidiness
problem.
"""

from __future__ import annotations

import argparse
import fnmatch
import getpass
import glob
import io
import os
import platform
import re
import socket
import sys
import tarfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Documents where third-party attribution legitimately carries e-mail addresses
# and personal names. Nothing else is allowed anywhere.
EMAIL_ALLOWED = [
    "LICENSE",
    "LICENSES/*",
    "NOTICE",
    "THIRD_PARTY_NOTICES.md",
    "THIRD_PARTY_LICENSES_RUST.md",
    "third_party/*",
    "*.dist-info/licenses/*",
    "*/LICENSE",
    "*/NOTICE",
]

SKIP_DIRS = {".git", "target", "__pycache__", ".venv", "venv", "node_modules", ".mypy_cache"}

# This file defines the patterns, so it contains one example of each by
# construction. Scanning it finds only itself. Anything else under tools/ is
# scanned normally.
SELF = os.path.join("tools", "privacy_scan.py")
BINARY_SUFFIXES = {".pyd", ".so", ".dll", ".exe", ".rlib", ".a", ".dylib", ".bin"}
ARCHIVE_SUFFIXES = {".zip", ".whl", ".crate", ".tar.gz", ".tgz", ".tar"}

EMAIL_RE = re.compile(rb"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")

PATTERNS = [
    ("win-path", re.compile(rb"(?i)[A-Za-z]:[\\/]+Users[\\/]+[A-Za-z0-9._-]+")),
    ("unix-path", re.compile(rb"(?:/home/|/Users/)[A-Za-z0-9._-]+")),
    ("unc", re.compile(rb"\\\\\\\\[A-Za-z0-9._-]+\\\\[A-Za-z0-9._$-]+")),
    ("env-path", re.compile(rb"(?i)%USERPROFILE%|\$HOME/|%HOMEPATH%")),
    ("toolchain-path", re.compile(rb"(?i)\.cargo[\\/]registry|\.rustup[\\/]toolchains")),
    ("email", EMAIL_RE),
]


def local_identity_pattern():
    """Names that identify *this* machine and account, resolved at run time.

    Deliberately not written down: putting the account name in the scanner would
    be the same mistake the scanner exists to catch.
    """
    tokens = set()
    for value in (
        getpass.getuser(),
        socket.gethostname(),
        platform.node(),
        os.environ.get("USERNAME", ""),
        os.environ.get("USER", ""),
        os.environ.get("COMPUTERNAME", ""),
    ):
        value = (value or "").strip()
        # Very short or generic tokens produce nothing but false positives.
        if len(value) >= 4 and value.lower() not in {"user", "root", "admin", "runner"}:
            tokens.add(value)
            if "." in value:
                tokens.add(value.split(".")[0])
    if not tokens:
        return None
    joined = b"|".join(re.escape(t.encode("utf-8", "ignore")) for t in sorted(tokens))
    # Word-boundary anchored, because an account name is often a short ordinary
    # string: a four-letter one matches inside "furnishing", "publishing" and
    # "distinguishing", which appear in every MIT and GPL text in the tree. An
    # unanchored match reports the licence files as leaking the user name and
    # buries the real findings under dozens of false ones.
    return re.compile(b"(?i)(?<![A-Za-z0-9])(?:" + joined + b")(?![A-Za-z0-9])")


def decodings(data):
    """The blob as bytes, plus a UTF-16LE-decoded view re-encoded to UTF-8.

    A Windows toolchain writes paths as UTF-16; searching only the raw bytes
    misses every one of them and reports a clean result on a dirty file.
    """
    yield "utf-8", data
    if b"\x00" in data[:4096]:
        try:
            yield "utf-16le", data.decode("utf-16le", "ignore").encode("utf-8", "ignore")
        except Exception:  # noqa: BLE001 - a malformed blob is not a finding
            pass


def retained_attribution_addresses(root=ROOT):
    """Every e-mail address the tree's own licence and notice documents carry.

    `src/licenses.rs` embeds those documents with `include_str!` so that
    `--licenses` can print them, which puts each upstream author's address inside
    every compiled binary. Those addresses are exactly what Apache-2.0 4(c) and
    the GPL require to be retained, so deleting them to quieten this scanner
    would break the licence terms it exists to protect.

    Membership in this set is the test, not the file the hit was found in. A
    binary is allowed to carry an address only because a licence document in the
    same tree carries it; an address that is in the binary and in no licence
    document is still a finding, which is what keeps this from becoming a blanket
    permit for binaries.
    """
    found = set()
    for pattern in EMAIL_ALLOWED:
        for path in glob.glob(os.path.join(root, pattern.replace("/", os.sep))):
            if not os.path.isfile(path):
                continue
            try:
                with open(path, "rb") as fh:
                    data = fh.read()
            except OSError:
                continue
            for match in EMAIL_RE.finditer(data):
                found.add(match.group(0).decode("utf-8", "replace").lower())
    return found


def allowed(name, cls, text="", retained=frozenset()):
    if cls != "email":
        return False
    posix = name.replace("\\", "/")
    if any(fnmatch.fnmatch(posix, pat) for pat in EMAIL_ALLOWED):
        return True
    return text.lower() in retained


def scan_blob(name, data, identity, findings, retained=frozenset()):
    if name.replace("\\", "/").endswith(SELF.replace("\\", "/")):
        return
    patterns = list(PATTERNS)
    if identity is not None:
        patterns.append(("local-identity", identity))
    for encoding, blob in decodings(data):
        for cls, pattern in patterns:
            for match in pattern.finditer(blob):
                text = match.group(0).decode("utf-8", "replace")
                if allowed(name, cls, text, retained):
                    continue
                findings.append((name, cls, encoding, text))
                if len(findings) > 500:
                    return


def scan_archive(path, identity, findings, retained=frozenset()):
    lower = path.lower()
    if lower.endswith((".zip", ".whl")):
        with zipfile.ZipFile(path) as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                scan_blob("{}!{}".format(os.path.basename(path), info.filename),
                          zf.read(info), identity, findings, retained)
    elif lower.endswith((".crate", ".tar.gz", ".tgz", ".tar")):
        mode = "r:gz" if lower.endswith((".crate", ".tar.gz", ".tgz")) else "r:"
        with tarfile.open(path, mode) as tf:
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                fh = tf.extractfile(member)
                if fh is None:
                    continue
                scan_blob("{}!{}".format(os.path.basename(path), member.name),
                          fh.read(), identity, findings, retained)
    else:
        with open(path, "rb") as fh:
            scan_blob(os.path.basename(path), fh.read(), identity, findings, retained)


def walk_tree(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            yield os.path.relpath(path, root), path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", metavar="DIR",
                    help="scan a source tree (skips .git, target, __pycache__)")
    ap.add_argument("--artifacts", nargs="*", metavar="FILE", default=[],
                    help="scan built artifacts: wheels, crates, binaries")
    args = ap.parse_args()
    if not args.tree and not args.artifacts:
        ap.error("give --tree and/or --artifacts")

    identity = local_identity_pattern()
    retained = retained_attribution_addresses()
    findings = []
    scanned = 0

    if args.tree:
        for name, path in walk_tree(args.tree):
            suffix = os.path.splitext(name)[1].lower()
            try:
                if suffix in ARCHIVE_SUFFIXES or name.lower().endswith(".tar.gz"):
                    scan_archive(path, identity, findings, retained)
                else:
                    with open(path, "rb") as fh:
                        scan_blob(name, fh.read(), identity, findings, retained)
                scanned += 1
            except OSError as exc:
                print("  could not read {}: {}".format(name, exc), file=sys.stderr)

    for pattern in args.artifacts:
        for path in sorted(glob.glob(pattern)):
            try:
                scan_archive(path, identity, findings, retained)
                scanned += 1
            except Exception as exc:  # noqa: BLE001 - report, do not mask
                print("  could not read {}: {}".format(path, exc), file=sys.stderr)
                findings.append((path, "unreadable", "-", str(exc)))

    if identity is None:
        print("note: no local identity tokens could be resolved on this machine; "
              "the local-identity class was not checked")

    if not findings:
        print("privacy scan: {} files, 0 findings".format(scanned))
        return 0

    print("privacy scan: {} files, {} finding(s)".format(scanned, len(findings)))
    by_class = {}
    for name, cls, encoding, text in findings:
        by_class.setdefault(cls, []).append((name, encoding, text))
    for cls in sorted(by_class):
        hits = by_class[cls]
        print("\n{} ({} hit(s)):".format(cls, len(hits)))
        seen = set()
        for name, encoding, text in hits:
            key = (name, text)
            if key in seen:
                continue
            seen.add(key)
            print("  {} [{}]  {}".format(name, encoding, text))
            if len(seen) >= 40:
                print("  ... and {} more in this class".format(len(hits) - len(seen)))
                break
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
