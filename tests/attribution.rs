// SPDX-License-Identifier: GPL-3.0-or-later
//! Pin the licensing and attribution record so it fails on drift.
//!
//! Licence compliance decays silently. A file gets embedded and nobody adds it
//! to the notices; a dependency is added and its licence is never recorded; a
//! well-meaning commit creates a NOTICE file upstream never wrote. None of that
//! breaks a build or a physics test, so nothing catches it.
//!
//! These tests catch it. Each one pins a specific obligation and names which
//! clause it answers, so a failure says what is now wrong rather than just that
//! a string changed.
//!
//! The v0.2.4 state these were written against had four real defects, and every
//! one of them is a test below:
//!
//! * `THIRD_PARTY_NOTICES.md` claimed only parameter data was taken, while a
//!   dozen modules are line-level ports of Apache-2.0 Fortran;
//! * no upstream copyright notice was retained anywhere (Apache-2.0 4(c));
//! * no modification notice or date (Apache-2.0 4(b), GPL-3.0 5(a));
//! * `src/data/element_data.csv` was embedded but listed nowhere.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(relative: &str) -> String {
    let path = root().join(relative);
    fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{} is required by the licence terms: {e}", path.display()))
}

fn sha256_hex(bytes: &[u8]) -> String {
    // A tiny SHA-256 so this test needs no dependency. The crate ships with
    // three dependencies and an attribution test is a poor reason for a fourth.
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let mut message = bytes.to_vec();
    let bit_len = (bytes.len() as u64) * 8;
    message.push(0x80);
    while message.len() % 64 != 56 {
        message.push(0);
    }
    message.extend_from_slice(&bit_len.to_be_bytes());
    for chunk in message.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[i * 4],
                chunk[i * 4 + 1],
                chunk[i * 4 + 2],
                chunk[i * 4 + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh) =
            (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        for (slot, value) in h.iter_mut().zip([a, b, c, d, e, f, g, hh]) {
            *slot = slot.wrapping_add(value);
        }
    }
    h.iter().map(|w| format!("{w:08x}")).collect()
}

/// Rust modules that port upstream code, and the upstream each derives from.
///
/// The second field is a substring that must appear in the file's `UPSTREAM:`
/// line **and** in `THIRD_PARTY_NOTICES.md`, so the two records cannot drift
/// apart.
const PORTED: &[(&str, &str)] = &[
    ("src/constants.rs", "conref_C.F90"),
    ("src/data_tables.rs", "parameters_C.F90"),
    ("src/integrals.rs", "mndod.F90"),
    ("src/integrals_d.rs", "mndod_C.F90"),
    ("src/onecenter.rs", "eiscor"),
    ("src/overlap.rs", "diat.F90"),
    ("src/params.rs", "calpar.F90"),
    ("src/repulsion.rs", "ccrep.F90"),
    ("src/rotations.rs", "rotmat"),
    ("src/scf.rs", "shell occupancies"),
    ("src/zindo.rs", "reimers_C.F90"),
    ("src/mindo3.rs", "block.f"),
    ("src/cndo_indo.rs", "Cndo2.cpp"),
    // Ported from gfn2-rs rather than from Fortran, so the marker is that
    // project's own file. Its numbers are Stewart's published tables, cited as
    // facts; its arrangement follows xtb's `slaterToGauss`, which is
    // LGPL-3.0-or-later. Both halves are in THIRD_PARTY_NOTICES.md, and this
    // row is what keeps the header on the file.
    ("src/sto.rs", "gfn2-rs"),
];

/// Data files embedded with `include_str!`. `element_data.csv` is the canary:
/// it was embedded and attributed nowhere through v0.2.4, which is the failure
/// mode these tests exist to prevent, so it is named explicitly rather than
/// left to a directory walk that could quietly find nothing.
const EMBEDDED_DATA: &[&str] = &[
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
];

#[test]
fn every_source_file_declares_its_licence() {
    let mut missing = Vec::new();
    for dir in ["src", "src/bin", "tests"] {
        let path = root().join(dir);
        let Ok(entries) = fs::read_dir(&path) else {
            continue;
        };
        for entry in entries.flatten() {
            let file = entry.path();
            if file.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let text = fs::read_to_string(&file).unwrap();
            let head: String = text.lines().take(5).collect::<Vec<_>>().join("\n");
            if !head.contains("SPDX-License-Identifier: GPL-3.0-or-later") {
                missing.push(file.display().to_string());
            }
        }
    }
    assert!(
        missing.is_empty(),
        "these files carry no SPDX header in their first five lines: {missing:#?}"
    );
}

#[test]
fn every_ported_module_carries_a_complete_provenance_header() {
    // Apache-2.0 4(b) and 4(c), and GPL-3.0 5(a): the attribution notice must be
    // retained and the fact and date of modification stated. A header that names
    // the upstream file but not the copyright holder satisfies neither.
    for (file, upstream_marker) in PORTED {
        let text = read(file);
        let head: String = text.lines().take(48).collect::<Vec<_>>().join("\n");
        assert!(
            head.contains("PROVENANCE:"),
            "{file} ports upstream code but has no PROVENANCE header"
        );
        assert!(
            head.contains("UPSTREAM:"),
            "{file} does not name the upstream file it derives from"
        );
        assert!(
            head.contains(upstream_marker),
            "{file} should name {upstream_marker:?} in its UPSTREAM line"
        );
        assert!(
            head.contains("MODIFIED for xndo-rs"),
            "{file} states no modification; GPL-3.0 5(a) requires it"
        );
        // An ISO date, so "when" is answerable.
        let has_date = head
            .split_whitespace()
            .any(|t| t.len() >= 10 && t.as_bytes()[4] == b'-' && t.as_bytes()[7] == b'-');
        assert!(has_date, "{file}'s MODIFIED line carries no ISO date");
        // Whose copyright, in the file itself.
        let attributed =
            head.contains("Copyright") || head.contains("public domain (no copyright asserted)");
        assert!(
            attributed,
            "{file} states no copyright holder; Apache-2.0 4(c) requires the \
             attribution notice to be retained, not merely referenced"
        );
    }
}

#[test]
fn the_notices_record_every_ported_module() {
    let notices = read("THIRD_PARTY_NOTICES.md");
    for (file, upstream_marker) in PORTED {
        assert!(
            notices.contains(file),
            "{file} is a port of upstream code but THIRD_PARTY_NOTICES.md does \
             not list it"
        );
        assert!(
            notices.contains(upstream_marker),
            "THIRD_PARTY_NOTICES.md does not say {file} derives from \
             {upstream_marker:?}"
        );
    }
    // The claim that was wrong through v0.2.4.
    assert!(
        notices.contains("line-level ports"),
        "THIRD_PARTY_NOTICES.md must say that code, not only parameter data, is \
         taken from OpenMOPAC; saying otherwise is what made v0.2.4's record \
         inaccurate"
    );
}

#[test]
fn every_embedded_data_file_is_attributed_and_hashed() {
    let notices = read("THIRD_PARTY_NOTICES.md");
    let manifest = read("src/data/legacy_parameter_manifest.sha256");
    for name in EMBEDDED_DATA {
        assert!(
            notices.contains(name),
            "src/data/{name} is embedded with include_str! but is absent from \
             THIRD_PARTY_NOTICES.md. This is exactly how element_data.csv went \
             unattributed through v0.2.4."
        );
        assert!(
            manifest.contains(name),
            "src/data/{name} is embedded but not hashed in \
             legacy_parameter_manifest.sha256, so nothing would notice it change"
        );
    }
}

#[test]
fn the_embedded_set_and_the_manifest_are_the_same_set() {
    // Catches a file added to one list and not the other, in either direction.
    let tables = read("src/data_tables.rs");
    let mut embedded = BTreeSet::new();
    for line in tables.lines() {
        if let Some(start) = line.find("include_str!(\"data/") {
            let rest = &line[start + "include_str!(\"data/".len()..];
            if let Some(end) = rest.find('"') {
                embedded.insert(rest[..end].to_string());
            }
        }
    }
    let manifest = read("src/data/legacy_parameter_manifest.sha256");
    let listed: BTreeSet<String> = manifest
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| l.split_whitespace().nth(1).map(str::to_string))
        .collect();
    // The manifest is itself embedded and cannot hash itself.
    let embedded: BTreeSet<String> = embedded
        .into_iter()
        .filter(|f| f != "legacy_parameter_manifest.sha256")
        .collect();
    assert_eq!(
        embedded, listed,
        "the set of include_str!-ed data files and the set in the manifest differ"
    );
}

#[test]
fn the_upstream_licence_copy_is_byte_identical() {
    // If this copy drifts from upstream's it is no longer the licence that was
    // granted, and Apache-2.0 4(a) is not satisfied by an approximation.
    let bytes = fs::read(root().join("third_party/mopac/LICENSE")).unwrap();
    assert_eq!(
        sha256_hex(&bytes),
        "6d1d968fb225eca367cb7f0b8831ab012a35d92b547e945e17ef8e7b05c3e5cc",
        "third_party/mopac/LICENSE no longer matches OpenMOPAC v23.2.5's copy"
    );
}

#[test]
fn the_retained_upstream_notices_are_present() {
    let notice = read("NOTICE");
    // Apache-2.0 4(c). This is the obligation v0.2.4 was failing outright.
    assert!(notice.contains("Copyright 2021"));
    assert!(notice.contains("Virginia Polytechnic Institute and State University"));
    assert!(notice.contains("Apache License, Version 2.0"));
    // MolDS, GPL-3.0-or-later.
    assert!(notice.contains("Mikiya Fujii"));
    // The INDO lineage upstream's own README records.
    assert!(notice.contains("Gieseking"));
    assert!(notice.contains("Reimers"));
}

#[test]
fn no_notice_file_is_invented_for_an_upstream_that_ships_none() {
    // Apache-2.0 4(d) obliges propagating a NOTICE only if the original work
    // has one. OpenMOPAC v23.2.5 does not. Putting a file at that path would
    // assert that it does, and would pass an obligation to carry a document
    // upstream never wrote on to every downstream recipient.
    let invented = root().join("third_party/mopac/NOTICE");
    assert!(
        !invented.exists(),
        "third_party/mopac/NOTICE must not exist: OpenMOPAC v23.2.5 ships no \
         NOTICE file, and creating one here would manufacture an obligation"
    );
    let provenance = read("third_party/mopac/PROVENANCE.md");
    assert!(
        provenance.contains("does not") || provenance.contains("ships no NOTICE"),
        "third_party/mopac/PROVENANCE.md must record that upstream ships no \
         NOTICE, so the absence is documented rather than looking like an \
         oversight"
    );
}

#[test]
fn every_referenced_licence_text_is_present_and_whole() {
    for id in [
        "Apache-2.0",
        "BSD-2-Clause",
        "BSD-3-Clause",
        "GPL-3.0-or-later",
        "MIT",
        "MPL-2.0",
        "Unicode-3.0",
        "Unlicense",
        "Zlib",
    ] {
        let text = read(&format!("LICENSES/{id}.txt"));
        assert!(
            text.len() > 500,
            "LICENSES/{id}.txt is {} bytes; a truncated licence is worse than a \
             missing one because it looks satisfied",
            text.len()
        );
    }
    // The project's own licence must be the same text in both places.
    let a = read("LICENSE")
        .replace("\r\n", "\n")
        .replace('\u{feff}', "");
    let b = read("LICENSES/GPL-3.0-or-later.txt").replace("\r\n", "\n");
    assert_eq!(
        a.trim(),
        b.trim(),
        "LICENSE and LICENSES/GPL-3.0-or-later.txt have diverged"
    );
}

#[test]
fn the_mpl_obligation_faer_hides_in_its_metadata_is_recorded() {
    // faer declares `license = "MIT"` and ships COPYING.EIGEN.MPL2 because its
    // divide-and-conquer SVD is ported from Eigen. A notices file generated from
    // SPDX metadata rather than from shipped files drops this, and faer is in
    // every artifact this project produces.
    let rust = read("THIRD_PARTY_LICENSES_RUST.md");
    assert!(rust.contains("COPYING.EIGEN.MPL2"));
    assert!(rust.contains("Mozilla Public License Version 2.0"));
    assert!(read("NOTICE").contains("COPYING.EIGEN.MPL2"));
    assert!(root().join("LICENSES/MPL-2.0.txt").exists());
}

#[test]
fn the_rust_notices_are_stamped_with_the_lockfile_they_describe() {
    let rust = read("THIRD_PARTY_LICENSES_RUST.md");
    let marker = "<!-- cargo-lock-sha256: ";
    let start = rust
        .find(marker)
        .expect("THIRD_PARTY_LICENSES_RUST.md carries no cargo-lock-sha256 marker");
    let recorded = &rust[start + marker.len()..start + marker.len() + 64];
    let actual = sha256_hex(&fs::read(root().join("Cargo.lock")).unwrap());
    assert_eq!(
        recorded, actual,
        "THIRD_PARTY_LICENSES_RUST.md was generated from a different Cargo.lock; \
         run `python tools/collect_rust_licenses.py`"
    );
}

#[test]
fn the_packaging_lists_carry_every_attribution_document() {
    // `include` in Cargo.toml is a whitelist and `license-files` in
    // pyproject.toml is what `pip show -f` and SBOM tools read. A document left
    // out of either is simply absent from that distribution, silently.
    let cargo = read("Cargo.toml");
    let pyproject = read("pyproject.toml");
    for required in [
        "NOTICE",
        "THIRD_PARTY_NOTICES.md",
        "THIRD_PARTY_LICENSES_RUST.md",
        "LICENSES/",
    ] {
        assert!(
            cargo.contains(required),
            "Cargo.toml `include` does not carry {required}, so it would be \
             absent from the published crate"
        );
        assert!(
            pyproject.contains(required),
            "pyproject.toml `license-files` does not carry {required}, so it \
             would be absent from the wheel"
        );
    }
    // Every licence text in the tree must be reachable by the glob.
    assert!(cargo.contains("/LICENSES/*.txt"));
    assert!(pyproject.contains("LICENSES/*.txt"));
}

#[test]
fn the_version_string_is_the_same_everywhere() {
    let cargo = read("Cargo.toml");
    let version = cargo
        .lines()
        .find(|l| l.starts_with("version"))
        .and_then(|l| l.split('=').nth(1))
        .map(|v| v.trim().trim_matches('"').to_string())
        .expect("Cargo.toml has no version");
    for (file, text) in [
        ("pyproject.toml", read("pyproject.toml")),
        (
            "python/xndo_rs/__init__.py",
            read("python/xndo_rs/__init__.py"),
        ),
        (
            "src/data/legacy_parameter_manifest.sha256",
            read("src/data/legacy_parameter_manifest.sha256"),
        ),
        (
            "src/data/legacy_parameter_catalog.csv",
            read("src/data/legacy_parameter_catalog.csv"),
        ),
    ] {
        assert!(
            text.contains(&version),
            "{file} does not mention version {version}; a version that is right \
             in one place and stale in another is how a release ships mislabelled"
        );
    }
}

#[test]
fn the_licenses_command_prints_the_documents() {
    // Both the library entry point and a real process, because the CLI wires
    // `--licenses` outside the option loop and that wiring is what would break.
    let documents = xndo_rs::third_party_licenses();
    assert!(documents.len() >= 12, "only {} documents", documents.len());
    assert!(documents.iter().any(|d| d.path == "NOTICE"));

    let output = std::process::Command::new(env!("CARGO_BIN_EXE_xndo_rs_cli"))
        .args(["licenses", "MPL"])
        .output()
        .expect("running the CLI");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("Mozilla Public License"),
        "`xndo_rs_cli licenses MPL` printed nothing useful"
    );

    // A bare invocation must print everything, not error.
    let all = std::process::Command::new(env!("CARGO_BIN_EXE_xndo_rs_cli"))
        .arg("licenses")
        .output()
        .expect("running the CLI");
    assert!(all.status.success());
    let all = String::from_utf8_lossy(&all.stdout);
    assert!(all.contains("Virginia Polytechnic Institute and State University"));
    assert!(all.len() > 100_000, "only {} bytes printed", all.len());
}

#[test]
fn the_privacy_scanner_exists_and_carries_no_identity_of_its_own() {
    // The scanner resolves the machine's identity at run time precisely so it
    // holds none itself. If someone hardcodes a user name into it, the tool
    // becomes the leak it was written to prevent.
    let path = root().join("tools/privacy_scan.py");
    assert!(path.exists(), "tools/privacy_scan.py is missing");
    let text = fs::read_to_string(&path).unwrap();
    assert!(text.contains("getpass.getuser()"));
    assert!(text.contains("socket.gethostname()"));
    // Assembled at run time rather than written out, because this file is
    // itself scanned: spelling the forbidden literals here would make the test
    // the thing the scanner reports. Found that out by running it.
    for (prefix, suffix) in [("C:\\Users", "\\s"), ("/home", "/s"), ("/Users", "/s")] {
        let literal = format!("{prefix}{suffix}");
        assert!(
            !text.contains(&literal),
            "tools/privacy_scan.py contains what looks like a real path \
             literal; the patterns must stay generic"
        );
    }
}

#[test]
fn upstream_attribution_documents_are_retained() {
    for (file, marker) in [
        ("third_party/mopac/AUTHORS.rst", "MOPAC Contributors"),
        ("third_party/mopac/CITATION.cff", "cff-version"),
        (
            "third_party/mopac/PROVENANCE.md",
            "052691223d19935a89f0fe18cd12301bd83e4201",
        ),
        ("third_party/molds/NOTICE", "Mikiya Fujii"),
        ("third_party/mopac7/NOTICE", "MOPAC7"),
    ] {
        let text = read(file);
        assert!(
            text.contains(marker),
            "{file} should contain {marker:?}; it is retained attribution and \
             must not be trimmed"
        );
    }
}

fn _unused(_: &Path) {}
