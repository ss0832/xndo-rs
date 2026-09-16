// SPDX-License-Identifier: GPL-3.0-or-later

use xndo_rs::data_tables::{
    CsvTable, LEGACY_PARAMETER_CATALOG_CSV, MINDO3_PAIR_PARAM_CSV, MINDO3_PARAM_CSV,
    MOLDS_CNDO2_INDO_PARAM_CSV, MOLDS_ZINDO_S_PARAM_CSV,
};

#[test]
fn bundled_legacy_parameter_tables_parse() {
    let cndo = CsvTable::parse(MOLDS_CNDO2_INDO_PARAM_CSV).unwrap();
    let zindo = CsvTable::parse(MOLDS_ZINDO_S_PARAM_CSV).unwrap();
    let m3 = CsvTable::parse(MINDO3_PARAM_CSV).unwrap();
    let m3p = CsvTable::parse(MINDO3_PAIR_PARAM_CSV).unwrap();
    let cat = CsvTable::parse(LEGACY_PARAMETER_CATALOG_CSV).unwrap();
    assert_eq!(cndo.rows.len(), 6);
    assert_eq!(zindo.rows.len(), 5);
    assert_eq!(m3.rows.len(), 10);
    assert_eq!(m3p.rows.len(), 40);
    // Seven since v0.3.0: `element_data.csv` was embedded but absent from the
    // catalog, the manifest and the notices, so nothing attributed or checked
    // it. See tests/attribution.rs.
    assert_eq!(cat.rows.len(), 7);
}

#[test]
fn molds_lithium_row_matches_gpl_source() {
    let table = CsvTable::parse(MOLDS_CNDO2_INDO_PARAM_CSV).unwrap();
    let z = table.col("z").unwrap();
    let imu_s = table.col("imu_amu_s_ev").unwrap();
    let g1 = table.col("indo_g1_native").unwrap();
    let row = table.rows.iter().find(|row| row[z] == "3").unwrap();
    assert_eq!(row[imu_s], "3.106");
    assert_eq!(row[g1], "0.092012");
}

#[test]
fn molds_sulfur_keeps_method_specific_d_and_zindo_values() {
    let cndo = CsvTable::parse(MOLDS_CNDO2_INDO_PARAM_CSV).unwrap();
    let z = cndo.col("z").unwrap();
    let d = cndo.col("has_d_cndo_indo").unwrap();
    let row = cndo.rows.iter().find(|r| r[z] == "16").unwrap();
    assert_eq!(row[d], "1");

    let zs = CsvTable::parse(MOLDS_ZINDO_S_PARAM_CSV).unwrap();
    let z = zs.col("z").unwrap();
    let beta = zs.col("beta_s_ev").unwrap();
    let row = zs.rows.iter().find(|r| r[z] == "16").unwrap();
    assert_eq!(row[beta], "-15.0");
}

#[test]
fn mindo3_chlorine_tail_matches_public_domain_mopac7() {
    let t = CsvTable::parse(MINDO3_PAIR_PARAM_CSV).unwrap();
    let z1 = t.col("z1").unwrap();
    let z2 = t.col("z2").unwrap();
    let beta = t.col("beta_ab").unwrap();
    let alpha = t.col("alpha_ab").unwrap();
    let p_cl = t
        .rows
        .iter()
        .find(|r| r[z1] == "15" && r[z2] == "17")
        .unwrap();
    let s_cl = t
        .rows
        .iter()
        .find(|r| r[z1] == "16" && r[z2] == "17")
        .unwrap();
    let cl_cl = t
        .rows
        .iter()
        .find(|r| r[z1] == "17" && r[z2] == "17")
        .unwrap();
    assert_eq!(
        (p_cl[beta].as_str(), p_cl[alpha].as_str()),
        ("0.277322", "1.543720")
    );
    assert_eq!(
        (s_cl[beta].as_str(), s_cl[alpha].as_str()),
        ("0.221764", "1.950318")
    );
    assert_eq!(
        (cl_cl[beta].as_str(), cl_cl[alpha].as_str()),
        ("0.258969", "1.792125")
    );
}
