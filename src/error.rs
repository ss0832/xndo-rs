// SPDX-License-Identifier: GPL-3.0-or-later

use std::fmt;

pub type Result<T> = std::result::Result<T, XndoError>;

/// Errors raised across the xndo-rs semiempirical pipeline.
///
/// The method is molecular NDDO: there is no periodic cell, and SCF
/// non-convergence is reported on the density residual.
#[derive(Debug)]
pub enum XndoError {
    Io(std::io::Error),
    Parse {
        line: usize,
        message: String,
    },
    InvalidInput(String),
    /// No parameter block exists for this atomic number in the selected method.
    MissingElement(u8),
    /// A code path was asked to handle a valence d shell that it does not support.
    /// The main MNDO pipeline supports d orbitals; this remains for deliberately
    /// s/p-only low-level helpers.
    UnsupportedDOrbitals(u8),
    /// A named per-element or derived parameter is absent.
    MissingParameter(String),
    LinearAlgebra(String),
    /// A requested calculation would exceed its configured soft workspace limit.
    ResourceLimit {
        operation: &'static str,
        required_mb: usize,
        limit_mb: usize,
    },
    /// The SCF loop hit `max_scf` without reaching the density/energy tolerance.
    ScfNotConverged {
        iterations: usize,
        error: f64,
    },
}

impl fmt::Display for XndoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "{err}"),
            Self::Parse { line, message } => write!(f, "parse error at line {line}: {message}"),
            Self::InvalidInput(msg) => write!(f, "{msg}"),
            Self::MissingElement(z) => write!(f, "missing semiempirical parameter block for Z={z}"),
            Self::UnsupportedDOrbitals(z) => {
                write!(f, "this calculation path does not support d orbitals for Z={z}")
            }
            Self::MissingParameter(key) => write!(f, "missing semiempirical parameter `{key}`"),
            Self::LinearAlgebra(msg) => write!(f, "linear algebra error: {msg}"),
            Self::ResourceLimit {
                operation,
                required_mb,
                limit_mb,
            } => write!(
                f,
                "{operation} requires approximately {required_mb} MiB, exceeding the configured {limit_mb} MiB workspace limit"
            ),
            Self::ScfNotConverged { iterations, error } => write!(
                f,
                "SCF did not converge after {iterations} iterations (error={error:.3e})"
            ),
        }
    }
}

impl std::error::Error for XndoError {}

impl From<std::io::Error> for XndoError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}
