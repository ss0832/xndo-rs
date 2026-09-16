// SPDX-License-Identifier: GPL-3.0-or-later

//! The licence and attribution documents, embedded in the library.
//!
//! Every obligation xndo-rs carries is discharged by a file in the source tree,
//! and every one of those files is embedded here with `include_str!`. That is
//! deliberate: a binary handed to someone without the source tree still carries
//! its own notices, and `tests/attribution.rs` can assert that the embedded set
//! and the tree agree rather than trusting that a packaging list was kept up to
//! date.
//!
//! Roughly 350 KB of the binary is these texts. That is the price of being able
//! to answer "what is in this?" from the artifact itself.

/// One embedded document: its path in the source tree and its contents.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LicenseDocument {
    /// Path relative to the repository root.
    pub path: &'static str,
    /// One line on what obligation this document answers.
    pub role: &'static str,
    pub text: &'static str,
}

/// The project's own licence.
pub const LICENSE: &str = include_str!("../LICENSE");

/// Retained upstream copyright notices (Apache-2.0 4(c), GPL-3.0 5(a)).
pub const NOTICE: &str = include_str!("../NOTICE");

/// The file-by-file record of what is taken from where.
pub const THIRD_PARTY_NOTICES: &str = include_str!("../THIRD_PARTY_NOTICES.md");

/// Rust crate licences, verbatim from what each crate ships.
pub const THIRD_PARTY_LICENSES_RUST: &str = include_str!("../THIRD_PARTY_LICENSES_RUST.md");

macro_rules! license_texts {
    ($(($id:literal, $file:literal)),* $(,)?) => {
        /// The licence texts referenced by this project or its dependencies.
        pub const LICENSE_TEXTS: &[(&str, &str)] = &[
            $(($id, include_str!(concat!("../LICENSES/", $file))),)*
        ];
    };
}

license_texts![
    ("Apache-2.0", "Apache-2.0.txt"),
    ("BSD-2-Clause", "BSD-2-Clause.txt"),
    ("BSD-3-Clause", "BSD-3-Clause.txt"),
    ("GPL-3.0-or-later", "GPL-3.0-or-later.txt"),
    ("LGPL-3.0-or-later", "LGPL-3.0-or-later.txt"),
    ("MIT", "MIT.txt"),
    ("MPL-2.0", "MPL-2.0.txt"),
    ("Unicode-3.0", "Unicode-3.0.txt"),
    ("Unlicense", "Unlicense.txt"),
    ("Zlib", "Zlib.txt"),
];

/// Every embedded attribution document, in the order `--licenses` prints them.
pub fn third_party_licenses() -> Vec<LicenseDocument> {
    let mut out = vec![
        LicenseDocument {
            path: "LICENSE",
            role: "the licence xndo-rs itself is distributed under",
            text: LICENSE,
        },
        LicenseDocument {
            path: "NOTICE",
            role: "retained upstream copyright notices",
            text: NOTICE,
        },
        LicenseDocument {
            path: "THIRD_PARTY_NOTICES.md",
            role: "what is taken from where, file by file",
            text: THIRD_PARTY_NOTICES,
        },
        LicenseDocument {
            path: "THIRD_PARTY_LICENSES_RUST.md",
            role: "Rust crate licences, verbatim from what each crate ships",
            text: THIRD_PARTY_LICENSES_RUST,
        },
    ];
    for (id, text) in LICENSE_TEXTS {
        out.push(LicenseDocument {
            path: "LICENSES/",
            role: id,
            text,
        });
    }
    out
}

/// Print the attribution documents.
///
/// With `filter`, prints only documents whose path or role contains it, so
/// `--licenses MPL` answers a question about one obligation without printing a
/// third of a megabyte.
pub fn print_licenses(filter: Option<&str>) {
    let documents = third_party_licenses();
    let matched: Vec<_> = documents
        .iter()
        .filter(|d| match filter {
            None => true,
            Some(f) => {
                let f = f.to_ascii_lowercase();
                d.path.to_ascii_lowercase().contains(&f) || d.role.to_ascii_lowercase().contains(&f)
            }
        })
        .collect();
    if matched.is_empty() {
        println!(
            "no attribution document matches {:?}. Available:",
            filter.unwrap_or("")
        );
        for d in &documents {
            println!("  {:<32} {}", d.path, d.role);
        }
        return;
    }
    for (index, d) in matched.iter().enumerate() {
        if index > 0 {
            println!();
        }
        println!("{}", "=".repeat(78));
        println!("{}  --  {}", d.path, d.role);
        println!("{}", "=".repeat(78));
        println!();
        println!("{}", d.text.trim_end());
    }
}

/// A one-line-per-document summary, for callers that want the list rather than
/// the texts.
pub fn license_summary() -> Vec<(&'static str, &'static str, usize)> {
    third_party_licenses()
        .into_iter()
        .map(|d| (d.path, d.role, d.text.len()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_document_is_embedded_and_non_trivial() {
        for d in third_party_licenses() {
            assert!(
                d.text.len() > 200,
                "{} ({}) is only {} bytes; an empty or truncated attribution \
                 document is worse than none, because it looks satisfied",
                d.path,
                d.role,
                d.text.len()
            );
        }
    }

    #[test]
    fn the_retained_upstream_copyright_is_present() {
        // Apache-2.0 4(c) is the obligation this project was failing before
        // v0.3.0: the notice has to be carried, not pointed at.
        assert!(NOTICE.contains("Copyright 2021"));
        assert!(NOTICE.contains("Virginia Polytechnic Institute and State University"));
        assert!(NOTICE.contains("Mikiya Fujii"));
    }

    #[test]
    fn the_mpl_obligation_faer_hides_in_its_metadata_is_recorded() {
        // faer declares MIT and ships COPYING.EIGEN.MPL2. If a future change
        // regenerates the notices with something that trusts the SPDX field,
        // this is what fails.
        assert!(LICENSE_TEXTS.iter().any(|(id, _)| *id == "MPL-2.0"));
        assert!(THIRD_PARTY_LICENSES_RUST.contains("COPYING.EIGEN.MPL2"));
        assert!(NOTICE.contains("COPYING.EIGEN.MPL2"));
    }

    #[test]
    fn filtering_finds_one_document_rather_than_all_of_them() {
        let all = third_party_licenses().len();
        let matched = third_party_licenses()
            .into_iter()
            .filter(|d| d.role.contains("MPL"))
            .count();
        assert_eq!(matched, 1);
        assert!(all > matched);
    }
}
