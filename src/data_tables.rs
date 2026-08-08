// SPDX-License-Identifier: GPL-3.0-or-later

//! Embedded semiempirical parameter/reference tables.
//!
//! Runtime MNDO, MNDO/d, and ZINDO/S tables come from OpenMOPAC v23.2.5
//! (Apache-2.0); MINDO/3 tables come from public-domain MOPAC7; reusable
//! CNDO2/INDO and ZINDO/S subsets come from MolDS 0.3.1
//! (GPL-3.0-or-later).
//! Every legacy CSV carries its own provenance header; see
//! `docs/parameter-provenance.md` and `THIRD_PARTY_NOTICES.md`.

/// MNDO per-element parameters from OpenMOPAC v23.2.5.
pub const MNDO_PARAM_CSV: &str = include_str!("data/mndo_parameters.csv");
/// MNDO pair-specific core-core parameters.
pub const MNDO_PAIR_CSV: &str = include_str!("data/mndo_pair_parameters.csv");
/// MNDO/d per-element parameters from OpenMOPAC v23.2.5.
pub const MNDOD_PARAM_CSV: &str = include_str!("data/mndod_parameters.csv");
/// MNDO/d pair-specific core-core parameters.
pub const MNDOD_PAIR_CSV: &str = include_str!("data/mndod_pair_parameters.csv");
/// Zerner INDO/S (ZINDO/S) spectroscopic parameters.
pub const ZINDO_S_PARAM_CSV: &str = include_str!("data/zindo_s_parameters.csv");
/// Reusable MolDS 0.3.1 atomic parameter subset used by its CNDO2/INDO paths
/// (H, Li, C, N, O, S; GPL-3.0-or-later upstream).
pub const MOLDS_CNDO2_INDO_PARAM_CSV: &str = include_str!("data/molds_cndo2_indo_parameters.csv");
/// Independent MolDS 0.3.1 ZINDO/S H/C/N/O/S parameter subset.
pub const MOLDS_ZINDO_S_PARAM_CSV: &str = include_str!("data/molds_zindo_s_parameters.csv");
/// Consolidated MINDO/3 per-element table cross-checked against public-domain MOPAC7.
pub const MINDO3_PARAM_CSV: &str = include_str!("data/mindo3_parameters.csv");
/// Public-domain MOPAC7 MINDO/3 atom-pair resonance/core-core table.
pub const MINDO3_PAIR_PARAM_CSV: &str = include_str!("data/mindo3_pair_parameters.csv");
/// Machine-readable inventory of all bundled legacy parameter datasets.
pub const LEGACY_PARAMETER_CATALOG_CSV: &str = include_str!("data/legacy_parameter_catalog.csv");
/// SHA-256 manifest freezing the exact xndo-rs 0.2.1 parameter snapshot.
pub const LEGACY_PARAMETER_MANIFEST_SHA256: &str =
    include_str!("data/legacy_parameter_manifest.sha256");

/// OpenMOPAC element reference data required by MNDO and MNDO/d
/// (occupancies, shell quantum numbers, atomic heats, and masses).
pub const ELEMENT_DATA_CSV: &str = include_str!("data/element_data.csv");

/// A parsed CSV: named columns and numeric rows (non-numeric cells become NaN,
/// the `sym` column is skipped by callers via name lookup on the raw fields).
pub struct CsvTable {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl CsvTable {
    /// Parse a comma-separated table, skipping `#` provenance/comment lines.
    pub fn parse(text: &str) -> Option<Self> {
        let mut lines = text.lines().filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        });
        let header = lines
            .next()?
            .split(',')
            .map(|s| s.trim().to_string())
            .collect::<Vec<_>>();
        let rows = lines
            .map(|l| l.split(',').map(|s| s.trim().to_string()).collect())
            .collect();
        Some(Self { header, rows })
    }

    pub fn col(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|c| c == name)
    }

    pub fn f64_at(&self, row: &[String], idx: usize) -> f64 {
        let value = row
            .get(idx)
            .unwrap_or_else(|| panic!("CSV row is missing column index {idx}"));
        if value.is_empty() {
            0.0
        } else {
            value
                .parse::<f64>()
                .unwrap_or_else(|_| panic!("invalid floating-point CSV value `{value}`"))
        }
    }
}

/// Per-element reference data parsed from `element_data.csv` (index = Z, 1..=107).
#[derive(Clone, Copy, Debug, Default)]
pub struct ElementData {
    /// Initial s/p/d shell occupancies (`ios`, `iop`, `iod`).
    pub occ_s: f64,
    pub occ_p: f64,
    pub occ_d: f64,
    /// Principal quantum numbers of the valence s/p/d shells (`npq`).
    pub npq_s: u8,
    pub npq_p: u8,
    pub npq_d: u8,
    /// Whether the one-center integrals derive directly from `Gss...Hsp` (true)
    /// or from `zsn/zpn/zdn` SlaterCondon parameters (false, transition metals).
    pub main_group: bool,
    /// d electrons assigned to the core in MOPAC's `Eisol` bookkeeping (`ndelec`).
    pub ndelec: i32,
    /// Experimental gas-phase atomic H_f (kcal/mol) (`eheat`).
    pub eheat_kcal: f64,
    /// Atomic mass (amu) (`ams`).
    pub mass: f64,
    /// Core charge = number of valence electrons (`tore`).
    pub tore: f64,
}

/// Parse `element_data.csv` into a Z-indexed table (index 0 unused).
pub fn element_data() -> Vec<ElementData> {
    let table = CsvTable::parse(ELEMENT_DATA_CSV).expect("embedded element_data.csv is valid");
    let col = |n: &str| {
        table
            .col(n)
            .unwrap_or_else(|| panic!("element_data.csv missing column {n}"))
    };
    let (c_z, c_ios, c_iop, c_iod) = (col("z"), col("ios"), col("iop"), col("iod"));
    let (c_ns, c_np, c_nd) = (col("npq_s"), col("npq_p"), col("npq_d"));
    let (c_mg, c_nde) = (col("main_group"), col("ndelec"));
    let (c_eh, c_mass, c_tore) = (col("eheat_kcal"), col("mass"), col("tore"));
    let mut out = vec![ElementData::default(); 108];
    for row in &table.rows {
        let z = table.f64_at(row, c_z) as usize;
        if z == 0 || z >= out.len() {
            continue;
        }
        out[z] = ElementData {
            occ_s: table.f64_at(row, c_ios),
            occ_p: table.f64_at(row, c_iop),
            occ_d: table.f64_at(row, c_iod),
            npq_s: table.f64_at(row, c_ns) as u8,
            npq_p: table.f64_at(row, c_np) as u8,
            npq_d: table.f64_at(row, c_nd) as u8,
            main_group: table.f64_at(row, c_mg) != 0.0,
            ndelec: table.f64_at(row, c_nde) as i32,
            eheat_kcal: table.f64_at(row, c_eh),
            mass: table.f64_at(row, c_mass),
            tore: table.f64_at(row, c_tore),
        };
    }
    out
}
