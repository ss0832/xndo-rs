// SPDX-License-Identifier: GPL-3.0-or-later

//! Dependency-light xndo-rs command-line frontend.

use std::env;
use std::path::Path;
use std::process::ExitCode;

use xndo_rs::constants::ANGSTROM_TO_BOHR;
use xndo_rs::{
    analytic_hessian, optimize, run_cndo_indo, run_gradient, run_hessian, run_mindo3,
    run_nddo_with_parameters, run_zindo_s, vibrational_analysis, vibrational_analysis_from_hessian,
    zindo_s_cis_gradients, zindo_s_cis_hessians, zindo_s_cis_spin, zindo_s_ucis,
    zindo_s_ucis_gradients, zindo_s_ucis_hessians, CisSpin, CndoIndoOptions, Method, Mindo3Options,
    Molecule, NddoOptions, NddoParameters, OptOptions, Reference, ZindoOptions, ZindoParameters,
};

fn cli_usage() -> &'static str {
    r#"Usage:
  xndo_rs_cli methods
  xndo_rs_cli <energy|gradient|charges|optimize|frequencies|hessian|uv-vis|spectrum|excited-properties|excited-gradient|excited-hessian> file.xyz [options]

Options:
  --method <name>                 Method (default: mndo)
  --charge <q>                    Molecular charge (default: 0)
  --multiplicity <m>              Spin multiplicity (default: 1)
  --reference <auto|rhf|uhf>      SCF reference for every native ground-state method
  --states <n>                    Number of ZINDO/S CIS roots per spin sector (default: 10)
  --state-type <singlet|triplet|both>
                                  ZINDO/S CIS spin sector (default: singlet)
  --no-diis                       Disable MNDO/MNDO-d DIIS
  --opt-max-iter <n>              Optimization iteration limit (default: 200)
  --opt-gtol <x>                  Max-gradient convergence in eV/Bohr (default: 1e-3)
  --opt-history <n>               L-BFGS history length (default: 8)

Executable method strings by API:
  energy, charges, single-point equivalent:
    CNDO/2   : cndo2 | cndo/2 | cndo
    INDO     : indo
    MNDO     : mndo
    MNDO/d   : mndod | mndo/d | mndo-d
    MINDO/3  : mindo3 | mindo/3 | mindo
    ZINDO/S  : zindo/s | zindos | zindo | indo/s | indos

  gradient, hessian, frequencies:
    CNDO/2   : cndo2 | cndo/2 | cndo
    INDO     : indo
    MNDO     : mndo
    MNDO/d   : mndod | mndo/d | mndo-d
    MINDO/3  : mindo3 | mindo/3 | mindo
    ZINDO/S  : zindo/s | zindos | zindo | indo/s | indos

  optimize:
    All six native ground-state methods listed above.

  uv-vis, spectrum, excited-properties, excited-gradient, excited-hessian:
    ZINDO/S  : zindo/s | zindos | zindo | indo/s | indos

RHF/UHF:
  Every native ground-state method accepts --reference auto|rhf|uhf.
  ZINDO/S excited-state commands use spin-adapted singlet/triplet RHF-CIS
  for RHF references and spin-orbital UCIS for UHF references.

Registered but parameter/engine-gated names:
  CNDO/1 INDO/1 INDO/2 ZINDO/1 ZINDO/2 MINDO/1 MINDO/2 SINDO1 MSINDO

CNDO/3, INDO/3 and ZINDO/3 are compatibility labels only and are rejected as
non-canonical methods. Run `xndo_rs_cli methods` for aliases and capabilities.
"#
}

#[derive(Debug)]
struct Cli {
    command: String,
    path: Option<String>,
    charge: f64,
    multiplicity: usize,
    use_diis: bool,
    method: Method,
    states: usize,
    reference: Reference,
    state_type: String,
    opt_max_iter: usize,
    opt_gtol: f64,
    opt_history: usize,
}

fn parse_reference(value: &str) -> Result<Reference, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "auto" | "" => Ok(Reference::Auto),
        "rhf" | "r" | "restricted" => Ok(Reference::Rhf),
        "uhf" | "u" | "unrestricted" => Ok(Reference::Uhf),
        _ => Err(format!(
            "reference must be auto, rhf, or uhf (got {value:?})"
        )),
    }
}

fn parse_args() -> Result<Cli, String> {
    let mut args = env::args().skip(1);
    let command = args.next().ok_or_else(|| cli_usage().to_owned())?;
    if command == "--help" || command == "-h" {
        return Err(String::new());
    }
    if command == "methods" {
        return Ok(Cli {
            command,
            path: None,
            charge: 0.0,
            multiplicity: 1,
            use_diis: true,
            method: Method::Mndo,
            states: 10,
            reference: Reference::Auto,
            state_type: "singlet".to_owned(),
            opt_max_iter: 200,
            opt_gtol: 1.0e-3,
            opt_history: 8,
        });
    }
    let path = args
        .next()
        .ok_or_else(|| format!("missing XYZ file\n\n{}", cli_usage()))?;
    let mut cli = Cli {
        command,
        path: Some(path),
        charge: 0.0,
        multiplicity: 1,
        use_diis: true,
        method: Method::Mndo,
        states: 10,
        reference: Reference::Auto,
        state_type: "singlet".to_owned(),
        opt_max_iter: 200,
        opt_gtol: 1.0e-3,
        opt_history: 8,
    };
    while let Some(flag) = args.next() {
        if flag == "--no-diis" {
            cli.use_diis = false;
            continue;
        }
        let value = args
            .next()
            .ok_or_else(|| format!("missing value for {flag}"))?;
        match flag.as_str() {
            "--charge" => {
                cli.charge = value
                    .parse()
                    .map_err(|_| format!("invalid charge: {value}"))?
            }
            "--multiplicity" => {
                cli.multiplicity = value
                    .parse()
                    .map_err(|_| format!("invalid multiplicity: {value}"))?;
                if cli.multiplicity == 0 {
                    return Err("multiplicity must be >= 1".into());
                }
            }
            "--method" => {
                cli.method = Method::parse(&value).ok_or_else(|| {
                    format!("unknown method {value:?}; run `xndo_rs_cli methods`")
                })?;
            }
            "--states" => {
                cli.states = value
                    .parse()
                    .map_err(|_| format!("invalid state count: {value}"))?
            }
            "--reference" => cli.reference = parse_reference(&value)?,
            "--opt-max-iter" => {
                cli.opt_max_iter = value
                    .parse()
                    .map_err(|_| format!("invalid optimization iteration limit: {value}"))?;
                if cli.opt_max_iter == 0 {
                    return Err("optimization iteration limit must be >= 1".into());
                }
            }
            "--opt-gtol" => {
                cli.opt_gtol = value
                    .parse()
                    .map_err(|_| format!("invalid optimization gradient tolerance: {value}"))?;
                if !cli.opt_gtol.is_finite() || cli.opt_gtol <= 0.0 {
                    return Err("optimization gradient tolerance must be finite and > 0".into());
                }
            }
            "--opt-history" => {
                cli.opt_history = value
                    .parse()
                    .map_err(|_| format!("invalid L-BFGS history length: {value}"))?;
                if cli.opt_history == 0 {
                    return Err("L-BFGS history length must be >= 1".into());
                }
            }
            "--state-type" => {
                let key = value.trim().to_ascii_lowercase();
                if !matches!(key.as_str(), "singlet" | "triplet" | "both") {
                    return Err(format!(
                        "state-type must be exactly singlet, triplet, or both (got {value:?})"
                    ));
                }
                cli.state_type = key;
            }
            _ => return Err(format!("unknown option: {flag}\n\n{}", cli_usage())),
        }
    }
    match cli.command.as_str() {
        "energy" | "gradient" | "charges" | "optimize" | "frequencies" | "hessian" | "uv-vis"
        | "spectrum" | "excited-properties" | "excited-gradient" | "excited-hessian" => {}
        _ => {
            return Err(format!(
                "unknown command: {}\n\n{}",
                cli.command,
                cli_usage()
            ))
        }
    }
    Ok(cli)
}

fn print_methods() {
    println!(
        "{:<10} {:<13} {:<13} {:<6} {:<6} {:<6} {:<8} {:<24} accepted strings",
        "method", "family", "status", "E", "grad", "Hess", "UHF", "APIs"
    );
    for method in Method::ALL {
        println!(
            "{:<10} {:<13} {:<13} {:<6} {:<6} {:<6} {:<8} {:<24} {}",
            method.as_str(),
            method.family(),
            method.status().as_str(),
            method.supports_energy(),
            method.supports_gradient(),
            method.supports_hessian(),
            method.supports_uhf(),
            method.api_names().join(","),
            method.accepted_strings().join("|")
        );
        println!("  note: {}", method.implementation_note());
    }
}

fn write_xyz(path: &Path, molecule: &Molecule, method: Method) -> std::io::Result<()> {
    let mut out = format!(
        "{}\nxndo-rs {} optimized geometry; coordinates in Angstrom\n",
        molecule.len(),
        method
    );
    for atom in &molecule.atoms {
        let symbol = xndo_rs::z_to_symbol(atom.z).unwrap_or("X");
        let pos = atom.position / ANGSTROM_TO_BOHR;
        out.push_str(&format!(
            "{symbol:2} {:+.10} {:+.10} {:+.10}\n",
            pos.x, pos.y, pos.z
        ));
    }
    std::fs::write(path, out)
}

fn nddo_options(cli: &Cli) -> NddoOptions {
    NddoOptions {
        charge: cli.charge,
        multiplicity: cli.multiplicity,
        reference: cli.reference,
        use_diis: cli.use_diis,
        ..NddoOptions::default()
    }
}

fn requested_spins(state_type: &str) -> Result<Vec<CisSpin>, String> {
    match state_type {
        "singlet" => Ok(vec![CisSpin::Singlet]),
        "triplet" => Ok(vec![CisSpin::Triplet]),
        "both" => Ok(vec![CisSpin::Singlet, CisSpin::Triplet]),
        other => Err(format!(
            "state-type must be exactly singlet, triplet, or both (got {other:?})"
        )),
    }
}

fn print_spectrum(
    molecule: &Molecule,
    params: &ZindoParameters,
    opts: &ZindoOptions,
    state_type: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let spectra = if opts.reference == Reference::Uhf || opts.multiplicity != 1 {
        vec![zindo_s_ucis(molecule, params, opts)?]
    } else {
        requested_spins(state_type)?
            .into_iter()
            .map(|spin| zindo_s_cis_spin(molecule, params, opts, spin))
            .collect::<xndo_rs::Result<Vec<_>>>()?
    };
    for spec in spectra {
        let spin = spec.spin;
        println!("# ZINDO/S {} CIS", spin.as_str());
        println!(
            "# ground dipole / D = {:+.8} {:+.8} {:+.8}",
            spec.ground_dipole_debye[0], spec.ground_dipole_debye[1], spec.ground_dipole_debye[2]
        );
        println!("# state spin E/eV cm^-1 lambda/nm f |mu_tr|/au |mu_exc|/D CT/A dominant");
        for (idx, st) in spec.states.iter().enumerate() {
            let dom = st
                .dominant
                .iter()
                .map(|c| {
                    format!(
                        "{}->{}:{:+.4}",
                        c.occupied, c.virtual_orbital, c.coefficient
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            println!(
                "{:>5} {:>7} {:>11.6} {:>12.2} {:>11.3} {:>11.7} {:>11.6} {:>11.6} {:>10.5} {}",
                idx + 1,
                st.spin.as_str(),
                st.energy_ev,
                st.energy_cm1,
                st.wavelength_nm.unwrap_or(f64::NAN),
                st.oscillator_strength,
                st.transition_dipole_magnitude_au,
                st.permanent_dipole_magnitude_debye,
                st.charge_transfer_distance_angstrom,
                dom,
            );
            println!(
                "#   transition dipole / D = {:+.8} {:+.8} {:+.8}",
                st.transition_dipole_debye[0],
                st.transition_dipole_debye[1],
                st.transition_dipole_debye[2]
            );
            println!(
                "#   spin multiplicity = {}; <S^2> = {:.6}; absolute state E / eV = {:.10}",
                st.spin_multiplicity, st.s2_expectation, st.state_total_energy_ev
            );
            println!(
                "#   excited permanent dipole / D = {:+.8} {:+.8} {:+.8}",
                st.permanent_dipole_debye[0],
                st.permanent_dipole_debye[1],
                st.permanent_dipole_debye[2]
            );
            println!(
                "#   difference dipole (exc-ground) / D = {:+.8} {:+.8} {:+.8}; |dmu| = {:.8}",
                st.difference_dipole_debye[0],
                st.difference_dipole_debye[1],
                st.difference_dipole_debye[2],
                st.difference_dipole_magnitude_debye
            );
            println!("#   hole centroid / A = {:+.8} {:+.8} {:+.8}; electron centroid / A = {:+.8} {:+.8} {:+.8}",
                st.hole_centroid_angstrom[0], st.hole_centroid_angstrom[1], st.hole_centroid_angstrom[2],
                st.electron_centroid_angstrom[0], st.electron_centroid_angstrom[1], st.electron_centroid_angstrom[2]);
            println!("#   state charges = {:?}", st.charges);
            println!("#   hole population = {:?}", st.hole_population);
            println!("#   electron population = {:?}", st.electron_population);
        }
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    if cli.command == "methods" {
        print_methods();
        return Ok(());
    }
    if !cli.method.supports_energy() {
        return Err(cli.method.execution_error().into());
    }
    let path = cli.path.as_deref().ok_or("missing XYZ file")?;
    let molecule = Molecule::from_xyz_file(path, cli.charge)?.with_multiplicity(cli.multiplicity);

    if cli.command == "optimize" {
        let result = optimize(
            &molecule,
            cli.method,
            &nddo_options(&cli),
            &OptOptions {
                max_iter: cli.opt_max_iter,
                gtol: cli.opt_gtol,
                history: cli.opt_history,
            },
        )?;
        let input = Path::new(path);
        let stem = input
            .file_stem()
            .and_then(|value| value.to_str())
            .unwrap_or("optimized");
        let output = input.with_file_name(format!(
            "{stem}.{}.opt.xyz",
            cli.method.as_str().to_ascii_lowercase().replace('/', "d")
        ));
        write_xyz(&output, &result.molecule, cli.method)?;
        println!(
            "Converged: {} after {} steps",
            result.converged, result.iterations
        );
        println!("Final energy: {:.12} eV", result.energy_ev);
        println!("Optimized geometry: {}", output.display());
        return Ok(());
    }

    if cli.command == "gradient" {
        let result = run_gradient(&molecule, cli.method, &nddo_options(&cli))?;
        println!("Method:                {}", cli.method);
        println!(
            "Reference:             {}",
            if result.unrestricted { "UHF" } else { "RHF" }
        );
        println!("Energy:                {:.12} eV", result.energy_ev);
        println!("# atom       dE/dx              dE/dy              dE/dz       (eV/Bohr)");
        for (i, g) in result.gradient.iter().enumerate() {
            println!(
                "{:>5} {:>18.10e} {:>18.10e} {:>18.10e}",
                i + 1,
                g.x,
                g.y,
                g.z
            );
        }
        return Ok(());
    }
    if cli.command == "hessian" || cli.command == "frequencies" {
        let result = run_hessian(&molecule, cli.method, &nddo_options(&cli))?;
        if cli.command == "hessian" {
            println!("Method:                {}", cli.method);
            println!("# Cartesian Hessian (eV/Bohr^2), row-major");
            for i in 0..result.hessian.rows {
                for j in 0..result.hessian.cols {
                    print!(" {:>18.10e}", result.hessian[(i, j)]);
                }
                println!();
            }
        } else {
            let vib = vibrational_analysis_from_hessian(&molecule, result.hessian)?;
            println!("Method:                {}", cli.method);
            println!("Frequencies (cm^-1):   {:?}", vib.frequencies_cm);
        }
        return Ok(());
    }

    match cli.method {
        Method::ZindoS => {
            let params = ZindoParameters::cached()?;
            let zopts = ZindoOptions {
                charge: cli.charge,
                multiplicity: cli.multiplicity,
                reference: cli.reference,
                n_states: cli.states,
                ..ZindoOptions::default()
            };
            match cli.command.as_str() {
                "energy" | "charges" => {
                    let r = run_zindo_s(&molecule, params, &zopts)?;
                    println!("Method:                ZINDO/S");
                    println!(
                        "Reference:             {}",
                        if r.unrestricted { "UHF" } else { "RHF" }
                    );
                    println!("Total energy:          {:.12} eV", r.total_ev);
                    println!("Electronic energy:     {:.12} eV", r.electronic_ev);
                    println!("Core-core energy:      {:.12} eV", r.core_ev);
                    println!("SCF iterations:        {}", r.iterations);
                    if cli.command == "charges" {
                        println!("# atom  element    INDO population charge (e)");
                        for (i, (atom, q)) in molecule.atoms.iter().zip(&r.charges).enumerate() {
                            println!(
                                "{:>5} {:>4} {:>20.10}",
                                i + 1,
                                xndo_rs::z_to_symbol(atom.z).unwrap_or("X"),
                                q
                            );
                        }
                    }
                }
                "uv-vis" | "spectrum" | "excited-properties" => {
                    print_spectrum(&molecule, params, &zopts, &cli.state_type)?;
                }
                "excited-gradient" => {
                    let results = if cli.reference == Reference::Uhf || cli.multiplicity != 1 {
                        vec![zindo_s_ucis_gradients(&molecule, params, &zopts)?]
                    } else {
                        requested_spins(&cli.state_type)?
                            .into_iter()
                            .map(|spin| zindo_s_cis_gradients(&molecule, params, &zopts, spin))
                            .collect::<xndo_rs::Result<Vec<_>>>()?
                    };
                    for result in results {
                        let spin = result.spin;
                        for state in result.states {
                            println!(
                                "# {} root {} state gradient (eV/Bohr)",
                                spin.as_str(),
                                state.root
                            );
                            for (atom, value) in state.state_gradient.iter().enumerate() {
                                println!(
                                    "{:>5} {:>18.10e} {:>18.10e} {:>18.10e}",
                                    atom + 1,
                                    value.x,
                                    value.y,
                                    value.z
                                );
                            }
                        }
                    }
                }
                "excited-hessian" => {
                    let results = if cli.reference == Reference::Uhf || cli.multiplicity != 1 {
                        vec![zindo_s_ucis_hessians(&molecule, params, &zopts)?]
                    } else {
                        requested_spins(&cli.state_type)?
                            .into_iter()
                            .map(|spin| zindo_s_cis_hessians(&molecule, params, &zopts, spin))
                            .collect::<xndo_rs::Result<Vec<_>>>()?
                    };
                    for result in results {
                        let spin = result.spin;
                        for state in result.states {
                            println!(
                                "# {} root {} state Hessian (eV/Bohr^2)",
                                spin.as_str(),
                                state.root
                            );
                            for i in 0..state.state_hessian.rows {
                                for j in 0..state.state_hessian.cols {
                                    print!(" {:>18.10e}", state.state_hessian[(i, j)]);
                                }
                                println!();
                            }
                        }
                    }
                }
                other => return Err(format!("ZINDO/S does not support `{other}`").into()),
            }
        }
        Method::Cndo2 | Method::Indo => match cli.command.as_str() {
            "energy" | "charges" => {
                let opts = CndoIndoOptions {
                    charge: cli.charge,
                    multiplicity: cli.multiplicity,
                    reference: cli.reference,
                    ..CndoIndoOptions::default()
                };
                let r = run_cndo_indo(&molecule, cli.method, &opts)?;
                println!("Method:                {}", cli.method);
                println!(
                    "Reference:             {}",
                    if r.unrestricted { "UHF" } else { "RHF" }
                );
                println!("N(alpha)/N(beta):      {}/{}", r.n_alpha, r.n_beta);
                println!("Total energy:          {:.12} eV", r.total_ev);
                println!("Electronic energy:     {:.12} eV", r.electronic_ev);
                println!("Core-core energy:      {:.12} eV", r.core_ev);
                println!("SCF iterations:        {}", r.iterations);
                println!("MO energies (eV):      {:?}", r.mo_energies_ev);
                if let Some(beta) = &r.mo_energies_beta_ev {
                    println!("Beta MO energies (eV): {:?}", beta);
                }
                if cli.command == "charges" {
                    println!("# atom  element    NDO population charge (e)");
                    for (i, (atom, q)) in molecule.atoms.iter().zip(&r.charges).enumerate() {
                        println!(
                            "{:>5} {:>4} {:>20.10}",
                            i + 1,
                            xndo_rs::z_to_symbol(atom.z).unwrap_or("X"),
                            q
                        );
                    }
                }
            }
            other => {
                return Err(format!(
                    "{} native path currently supports energy/charges only, not `{other}`",
                    cli.method
                )
                .into())
            }
        },
        Method::Mindo3 => match cli.command.as_str() {
            "energy" | "charges" => {
                let opts = Mindo3Options {
                    charge: cli.charge,
                    multiplicity: cli.multiplicity,
                    reference: cli.reference,
                    ..Mindo3Options::default()
                };
                let r = run_mindo3(&molecule, &opts)?;
                println!("Method:                MINDO/3");
                println!(
                    "Reference:             {}",
                    if r.unrestricted { "UHF" } else { "RHF" }
                );
                println!("N(alpha)/N(beta):      {}/{}", r.n_alpha, r.n_beta);
                println!("Total energy:          {:.12} eV", r.total_ev);
                println!("Electronic energy:     {:.12} eV", r.electronic_ev);
                println!("Core-core energy:      {:.12} eV", r.core_ev);
                println!(
                    "Heat of formation:     {:.8} kcal/mol",
                    r.heat_of_formation_kcal
                );
                println!("SCF iterations:        {}", r.iterations);
                if cli.command == "charges" {
                    println!("# atom  element    population charge (e)");
                    for (i, (atom, q)) in molecule.atoms.iter().zip(&r.charges).enumerate() {
                        println!(
                            "{:>5} {:>4} {:>20.10}",
                            i + 1,
                            xndo_rs::z_to_symbol(atom.z).unwrap_or("X"),
                            q
                        );
                    }
                }
            }
            other => {
                return Err(format!(
                    "MINDO/3 native path currently supports energy/charges only, not `{other}`"
                )
                .into())
            }
        },
        Method::Mndo | Method::MndoD => {
            if matches!(
                cli.command.as_str(),
                "uv-vis" | "spectrum" | "excited-properties"
            ) {
                return Err("excited-state properties currently require --method zindo/s".into());
            }
            let params = NddoParameters::for_method(cli.method)?;
            let options = nddo_options(&cli);
            match cli.command.as_str() {
                "energy" => {
                    let r = run_nddo_with_parameters(&molecule, &params, &options)?;
                    println!("Method:                {}", cli.method);
                    println!(
                        "Reference:             {}",
                        if r.unrestricted { "UHF" } else { "RHF" }
                    );
                    println!("Total energy:          {:.12} eV", r.total_ev);
                    println!("Electronic energy:     {:.12} eV", r.electronic_ev);
                    println!("Core-core energy:      {:.12} eV", r.core_ev);
                    println!(
                        "Heat of formation:     {:.8} kcal/mol",
                        r.heat_of_formation_kcal
                    );
                    println!("SCF iterations:        {}", r.iterations);
                    println!("MO energies (eV):      {:?}", r.mo_energies);
                }
                "charges" => {
                    let r = run_nddo_with_parameters(&molecule, &params, &options)?;
                    println!("Reference: {}", if r.unrestricted { "UHF" } else { "RHF" });
                    println!("# atom  element    Mulliken charge (e)");
                    for (i, (atom, q)) in molecule.atoms.iter().zip(&r.charges).enumerate() {
                        println!(
                            "{:>5} {:>4} {:>20.10}",
                            i + 1,
                            xndo_rs::z_to_symbol(atom.z).unwrap_or("X"),
                            q
                        );
                    }
                    println!(
                        "Dipole: {:.8} {:.8} {:.8} D; |mu|={:.8} D",
                        r.dipole_debye.x, r.dipole_debye.y, r.dipole_debye.z, r.dipole_magnitude
                    );
                }
                "frequencies" => {
                    let r = vibrational_analysis(&molecule, &params, &options, 1.0e-3)?;
                    for (i, f) in r.frequencies_cm.iter().enumerate() {
                        println!("{:>5} {:>20.8}", i + 1, f);
                    }
                }
                "hessian" => {
                    let h = analytic_hessian(&molecule, &params, &options, 1.0e-3)?;
                    for i in 0..h.rows {
                        for j in 0..h.cols {
                            if j > 0 {
                                print!(" ");
                            }
                            print!("{:.10e}", h[(i, j)]);
                        }
                        println!();
                    }
                }
                _ => unreachable!(),
            }
        }
        other => return Err(other.execution_error().into()),
    }
    Ok(())
}

fn main() -> ExitCode {
    if matches!(
        env::args().nth(1).as_deref(),
        Some("--version") | Some("-V")
    ) {
        println!("xndo-rs {}", env!("CARGO_PKG_VERSION"));
        return ExitCode::SUCCESS;
    }
    match parse_args() {
        Ok(cli) => match run(cli) {
            Ok(()) => ExitCode::SUCCESS,
            Err(err) => {
                eprintln!("xndo-rs: {err}");
                ExitCode::FAILURE
            }
        },
        Err(message) if message.is_empty() => {
            print!("{}", cli_usage());
            ExitCode::SUCCESS
        }
        Err(message) => {
            eprintln!("xndo-rs: {message}");
            ExitCode::FAILURE
        }
    }
}
