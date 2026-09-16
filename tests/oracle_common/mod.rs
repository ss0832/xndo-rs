// SPDX-License-Identifier: GPL-3.0-or-later
//! Loader for the TSV oracle reference files in `tests/data/`.
//!
//! One column-name-keyed loader serves every reference file, so that a
//! multiplicity of files does not become a multiplicity of parsers. The format
//! mirrors the `#`-commented, header-carrying CSV convention the shipped
//! parameter tables already use (`src/data_tables.rs`), with tabs instead of
//! commas and an extra provenance block; a reviewer learns one format.
//!
//! Encodings, fixed project-wide:
//!
//! * a list cell is comma-separated;
//! * a geometry cell is `|`-separated, with `Sym x y z` inside, in Angstrom;
//! * an **empty cell means the oracle did not report that quantity**, and the
//!   comparison is skipped rather than run against a fabricated zero.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::PathBuf;

/// One cell of a reference row.
#[derive(Clone, Debug, PartialEq)]
pub enum Cell {
    /// The oracle did not report this quantity.
    Absent,
    Scalar(f64),
    List(Vec<f64>),
    Text(String),
}

impl Cell {
    fn parse(raw: &str) -> Self {
        let raw = raw.trim();
        if raw.is_empty() {
            return Cell::Absent;
        }
        if raw.contains(',') {
            let parts: Vec<&str> = raw.split(',').collect();
            let mut values = Vec::with_capacity(parts.len());
            for part in &parts {
                match part.trim().parse::<f64>() {
                    Ok(v) => values.push(v),
                    Err(_) => return Cell::Text(raw.to_string()),
                }
            }
            return Cell::List(values);
        }
        match raw.parse::<f64>() {
            Ok(v) => Cell::Scalar(v),
            Err(_) => Cell::Text(raw.to_string()),
        }
    }
}

/// Where a reference file came from, asserted by the harness so that swapping
/// the oracle for a later release fails loudly instead of silently.
#[derive(Clone, Debug, Default)]
pub struct Provenance {
    pub oracle_program: String,
    pub oracle_version: String,
    pub oracle_keywords: String,
    pub generator: String,
    pub generator_sha256: String,
    pub generated_utc: String,
}

#[derive(Clone, Debug)]
pub struct Row {
    pub name: String,
    pub charge: f64,
    pub multiplicity: usize,
    /// `(symbol, x, y, z)` in Angstrom, exactly as sent to the oracle.
    pub atoms: Vec<(String, f64, f64, f64)>,
    pub why: String,
    cells: BTreeMap<String, Cell>,
}

impl Row {
    /// XYZ text ready for `Molecule::from_xyz_str`.
    pub fn xyz(&self) -> String {
        let mut out = format!("{}\n{}\n", self.atoms.len(), self.name);
        for (sym, x, y, z) in &self.atoms {
            out.push_str(&format!("{sym} {x:.10} {y:.10} {z:.10}\n"));
        }
        out
    }

    pub fn elements(&self) -> Vec<String> {
        let mut seen: Vec<String> = self.atoms.iter().map(|a| a.0.clone()).collect();
        seen.sort();
        seen.dedup();
        seen
    }

    /// `None` when the oracle did not report the quantity.
    pub fn scalar(&self, column: &str) -> Option<f64> {
        match self.cells.get(column) {
            Some(Cell::Scalar(v)) => Some(*v),
            _ => None,
        }
    }

    /// Empty when the oracle did not report the quantity. A single-element list
    /// arrives as a `Scalar`, so accept both.
    pub fn list(&self, column: &str) -> Vec<f64> {
        match self.cells.get(column) {
            Some(Cell::List(v)) => v.clone(),
            Some(Cell::Scalar(v)) => vec![*v],
            _ => Vec::new(),
        }
    }

    pub fn has(&self, column: &str) -> bool {
        !matches!(
            self.cells.get(column),
            None | Some(Cell::Absent) | Some(Cell::Text(_))
        )
    }
}

#[derive(Clone, Debug)]
pub struct ReferenceFile {
    pub path: PathBuf,
    pub provenance: Provenance,
    pub rows: Vec<Row>,
}

fn parse_geometry(field: &str) -> Vec<(String, f64, f64, f64)> {
    field
        .split('|')
        .filter(|s| !s.trim().is_empty())
        .map(|atom| {
            let parts: Vec<&str> = atom.split_whitespace().collect();
            assert_eq!(
                parts.len(),
                4,
                "geometry atom {atom:?} does not have four fields"
            );
            (
                parts[0].to_string(),
                parts[1].parse().expect("geometry x"),
                parts[2].parse().expect("geometry y"),
                parts[3].parse().expect("geometry z"),
            )
        })
        .collect()
}

/// Load `tests/data/<name>`.
///
/// Panics with the regeneration command when the file is missing, because a
/// silently skipped oracle suite is worse than a failing one.
pub fn load(name: &str) -> ReferenceFile {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|err| {
        panic!(
            "cannot read oracle reference {}: {err}\n\
             regenerate it with:\n  \
             set MOPAC_EXE=<path to mopac.exe>\n  \
             python tools/oracle/build_reference_set.py --method all",
            path.display()
        )
    });

    let mut provenance = Provenance::default();
    let mut header: Vec<String> = Vec::new();
    let mut rows = Vec::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('#') {
            let mut it = rest.trim().splitn(2, '\t');
            let key = it.next().unwrap_or("").trim();
            let value = it.next().unwrap_or("").trim().to_string();
            match key {
                "oracle_program" => provenance.oracle_program = value,
                "oracle_version" => provenance.oracle_version = value,
                "oracle_keywords" => provenance.oracle_keywords = value,
                "generator" => provenance.generator = value,
                "generator_sha256" => provenance.generator_sha256 = value,
                "generated_utc" => provenance.generated_utc = value,
                _ => {}
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if header.is_empty() {
            header = fields.iter().map(|s| s.trim().to_string()).collect();
            continue;
        }
        assert_eq!(
            fields.len(),
            header.len(),
            "{}: row {:?} has {} fields, header has {}",
            path.display(),
            fields.first().unwrap_or(&""),
            fields.len(),
            header.len()
        );
        let mut cells = BTreeMap::new();
        for (column, raw) in header.iter().zip(&fields) {
            cells.insert(column.clone(), Cell::parse(raw));
        }
        let get_text = |column: &str| -> String {
            fields[header
                .iter()
                .position(|h| h == column)
                .unwrap_or_else(|| panic!("{}: missing required column {column}", path.display()))]
            .trim()
            .to_string()
        };
        rows.push(Row {
            name: get_text("name"),
            charge: get_text("charge").parse().expect("charge"),
            multiplicity: get_text("multiplicity").parse().expect("multiplicity"),
            atoms: parse_geometry(&get_text("geometry")),
            why: get_text("why"),
            cells,
        });
    }

    assert!(
        !rows.is_empty(),
        "{} contains no reference rows",
        path.display()
    );
    ReferenceFile {
        path,
        provenance,
        rows,
    }
}

/// Worst absolute deviation over a set of comparisons, with the row and column
/// that produced it, so a failure names the molecule rather than a number.
#[derive(Clone, Debug, Default)]
pub struct Worst {
    pub value: f64,
    pub label: String,
    pub compared: usize,
}

impl Worst {
    pub fn observe(&mut self, deviation: f64, label: impl Into<String>) {
        self.compared += 1;
        if deviation > self.value || self.label.is_empty() {
            self.value = deviation;
            self.label = label.into();
        }
    }
}
