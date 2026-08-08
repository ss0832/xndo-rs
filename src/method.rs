// SPDX-License-Identifier: GPL-3.0-or-later

//! Semiempirical Hamiltonians exposed by xndo-rs.
//!
//! xndo-rs deliberately distinguishes a method name from its implementation
//! status.  Several 1960s/1970s labels have been used inconsistently in old
//! program manuals.  Canonical names are therefore registered explicitly and
//! non-canonical compatibility labels (CNDO/3, INDO/3, ZINDO/3) fail loudly
//! instead of silently selecting a different Hamiltonian.

use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MethodStatus {
    /// Native numerical implementation is present in this source tree.
    Native,
    /// The canonical method is registered, with API/parameter hooks retained,
    /// but a provenance-complete built-in parameterization is not shipped yet.
    Registered,
    /// Requested compatibility spelling for which no canonical historical
    /// method definition was found.  Execution is intentionally rejected.
    NonCanonical,
}

impl MethodStatus {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Native => "native",
            Self::Registered => "registered",
            Self::NonCanonical => "noncanonical",
        }
    }
}

/// Semiempirical model selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Method {
    Cndo1,
    Cndo2,
    /// Compatibility label only; CNDO/3 is not a canonical Pople method.
    Cndo3,
    /// Unversioned ground-state INDO implementation label used by MolDS.
    Indo,
    Indo1,
    Indo2,
    /// Compatibility label only; INDO/3 is not a canonical Pople method.
    Indo3,
    Zindo1,
    Zindo2,
    /// Compatibility label only; ZINDO/3 is not a canonical Zerner method.
    Zindo3,
    /// Zerner INDO/S (ZINDO/S), fixed-geometry spectroscopy.
    ZindoS,
    Mindo1,
    Mindo2,
    Mindo3,
    Sindo,
    Msindo,
    /// Dewar-Thiel MNDO (NDDO).
    #[default]
    Mndo,
    /// Thiel-Voityuk MNDO/d (NDDO with explicit d functions).
    MndoD,
}

impl Method {
    pub const ALL: [Method; 18] = [
        Self::Cndo1,
        Self::Cndo2,
        Self::Cndo3,
        Self::Indo,
        Self::Indo1,
        Self::Indo2,
        Self::Indo3,
        Self::Zindo1,
        Self::Zindo2,
        Self::Zindo3,
        Self::ZindoS,
        Self::Mindo1,
        Self::Mindo2,
        Self::Mindo3,
        Self::Sindo,
        Self::Msindo,
        Self::Mndo,
        Self::MndoD,
    ];

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Cndo1 => "CNDO/1",
            Self::Cndo2 => "CNDO/2",
            Self::Cndo3 => "CNDO/3",
            Self::Indo => "INDO",
            Self::Indo1 => "INDO/1",
            Self::Indo2 => "INDO/2",
            Self::Indo3 => "INDO/3",
            Self::Zindo1 => "ZINDO/1",
            Self::Zindo2 => "ZINDO/2",
            Self::Zindo3 => "ZINDO/3",
            Self::ZindoS => "ZINDO/S",
            Self::Mindo1 => "MINDO/1",
            Self::Mindo2 => "MINDO/2",
            Self::Mindo3 => "MINDO/3",
            Self::Sindo => "SINDO1",
            Self::Msindo => "MSINDO",
            Self::Mndo => "MNDO",
            Self::MndoD => "MNDO/d",
        }
    }

    pub const fn family(self) -> &'static str {
        match self {
            Self::Cndo1 | Self::Cndo2 | Self::Cndo3 => "CNDO",
            Self::Indo | Self::Indo1 | Self::Indo2 | Self::Indo3 => "INDO",
            Self::Zindo1 | Self::Zindo2 | Self::Zindo3 | Self::ZindoS => "ZINDO/INDO",
            Self::Mindo1 | Self::Mindo2 | Self::Mindo3 => "MINDO",
            Self::Sindo => "SINDO",
            Self::Msindo => "MSINDO",
            Self::Mndo | Self::MndoD => "MNDO/NDDO",
        }
    }

    pub const fn status(self) -> MethodStatus {
        match self {
            Self::Cndo2 | Self::Indo | Self::Mndo | Self::MndoD | Self::Mindo3 | Self::ZindoS => {
                MethodStatus::Native
            }
            Self::Cndo3 | Self::Indo3 | Self::Zindo3 => MethodStatus::NonCanonical,
            _ => MethodStatus::Registered,
        }
    }

    pub const fn supports_energy(self) -> bool {
        matches!(
            self,
            Self::Cndo2 | Self::Indo | Self::Mndo | Self::MndoD | Self::Mindo3 | Self::ZindoS
        )
    }

    pub const fn supports_gradient(self) -> bool {
        matches!(
            self,
            Self::Cndo2 | Self::Indo | Self::Mndo | Self::MndoD | Self::Mindo3 | Self::ZindoS
        )
    }

    pub const fn supports_spectrum(self) -> bool {
        matches!(self, Self::ZindoS)
    }

    pub const fn supports_hessian(self) -> bool {
        matches!(
            self,
            Self::Cndo2 | Self::Indo | Self::Mndo | Self::MndoD | Self::Mindo3 | Self::ZindoS
        )
    }

    pub const fn supports_uhf(self) -> bool {
        matches!(
            self,
            Self::Cndo2 | Self::Indo | Self::Mndo | Self::MndoD | Self::Mindo3 | Self::ZindoS
        )
    }

    pub const fn supports_excited_properties(self) -> bool {
        matches!(self, Self::ZindoS)
    }

    /// Accepted user-facing method strings. Parsing is ASCII case-insensitive;
    /// spaces, `_`, and `-` are ignored by [`Method::parse`]. This list is
    /// intentionally explicit so Python/Rust/CLI callers can discover aliases
    /// without reading parser source.
    pub const fn accepted_strings(self) -> &'static [&'static str] {
        match self {
            Self::Cndo1 => &["cndo1", "cndo/1"],
            Self::Cndo2 => &["cndo2", "cndo/2", "cndo"],
            Self::Cndo3 => &["cndo3", "cndo/3"],
            Self::Indo => &["indo"],
            Self::Indo1 => &["indo1", "indo/1"],
            Self::Indo2 => &["indo2", "indo/2"],
            Self::Indo3 => &["indo3", "indo/3"],
            Self::Zindo1 => &["zindo1", "zindo/1"],
            Self::Zindo2 => &["zindo2", "zindo/2"],
            Self::Zindo3 => &["zindo3", "zindo/3"],
            Self::ZindoS => &["zindo/s", "zindos", "zindo", "indo/s", "indos"],
            Self::Mindo1 => &["mindo1", "mindo/1"],
            Self::Mindo2 => &["mindo2", "mindo/2"],
            Self::Mindo3 => &["mindo3", "mindo/3", "mindo"],
            Self::Sindo => &["sindo1", "sindo/1", "sindo"],
            Self::Msindo => &["msindo"],
            Self::Mndo => &["mndo"],
            Self::MndoD => &["mndod", "mndo/d", "mndo-d"],
        }
    }

    /// Public computational APIs that can execute this method in v0.2.4.
    pub const fn api_names(self) -> &'static [&'static str] {
        match self {
            Self::Cndo2 | Self::Indo | Self::Mindo3 => &[
                "single_point",
                "gradient",
                "forces",
                "optimize",
                "hessian",
                "frequencies",
            ],
            Self::Mndo | Self::MndoD => &[
                "single_point",
                "gradient",
                "forces",
                "optimize",
                "hessian",
                "frequencies",
            ],
            Self::ZindoS => &[
                "single_point",
                "gradient",
                "forces",
                "optimize",
                "hessian",
                "frequencies",
                "excited_states",
                "excited_properties",
                "uv_vis_spectrum",
                "excited_state_gradients",
                "excited_state_hessians",
            ],
            _ => &[],
        }
    }

    pub const fn implementation_note(self) -> &'static str {
        match self {
            Self::Mndo => "native NDDO SCF; OpenMOPAC MNDO parameter set",
            Self::MndoD => "native spd NDDO SCF; OpenMOPAC MNDOD parameter set",
            Self::Mindo3 => "native RHF/UHF MINDO/3 with analytic ground-state derivatives; canonical main-group parameter set",
            Self::ZindoS => "native RHF/UHF ground-state derivatives; spin-adapted RHF-CIS and spin-orbital UCIS properties and state derivatives for the s/p branch; transition-metal 9-AO branch is intentionally gated",
            Self::Cndo1 => "canonical method registered; awaiting provenance-complete CNDO/1 parameter/integral set",
            Self::Cndo2 => "native MolDS-compatible CNDO/2 RHF/UHF engine with analytic ground-state derivatives; bundled GPL-3.0-or-later H/Li/C/N/O/S donor parameters",
            Self::Indo => "native MolDS-compatible ground-state INDO RHF/UHF engine with analytic derivatives for bundled H/Li/C/N/O parameters; intentionally distinct from ZINDO/S (INDO/S)",
            Self::Indo1 => "canonical method registered; awaiting provenance-complete INDO/1 Slater-Condon parameter set",
            Self::Indo2 => "canonical method registered; awaiting provenance-complete INDO/2 Slater-Condon parameter set",
            Self::Zindo1 => "canonical Zerner ground-state model registered; built-in ZINDO/1 parameterization not yet shipped",
            Self::Zindo2 => "canonical Zerner model registered; built-in ZINDO/2 parameterization not yet shipped",
            Self::Mindo1 => "canonical MINDO/1 registered; complete validated pair parameterization not yet shipped",
            Self::Mindo2 => "canonical MINDO/2 registered; complete validated pair parameterization not yet shipped",
            Self::Sindo => "SINDO1 registered; orthogonalization/pseudopotential parameter set not redistributed",
            Self::Msindo => "MSINDO registered; no complete openly redistributable parameter and Hamiltonian set is bundled",
            Self::Cndo3 => "no canonical CNDO/3 definition identified; use CNDO/1 or CNDO/2",
            Self::Indo3 => "no canonical INDO/3 definition identified; use INDO/1, INDO/2 or INDO/S",
            Self::Zindo3 => "no canonical ZINDO/3 definition identified; use ZINDO/1, ZINDO/2 or ZINDO/S",
        }
    }

    pub fn execution_error(self) -> String {
        match self.status() {
            MethodStatus::Native => format!("{} is native and should be dispatched by its engine", self),
            MethodStatus::Registered => format!(
                "{} is registered in xndo-rs, but this release does not yet expose a complete executable numerical implementation: {}",
                self, self.implementation_note()
            ),
            MethodStatus::NonCanonical => format!(
                "{} is retained only as a compatibility/requested label and is not a canonical historical method: {}",
                self, self.implementation_note()
            ),
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        let key = s.trim().to_ascii_lowercase().replace(['_', ' ', '-'], "");
        match key.as_str() {
            "cndo1" | "cndo/1" => Some(Self::Cndo1),
            "cndo2" | "cndo/2" | "cndo" => Some(Self::Cndo2),
            "cndo3" | "cndo/3" => Some(Self::Cndo3),
            "indo" => Some(Self::Indo),
            "indo1" | "indo/1" => Some(Self::Indo1),
            "indo2" | "indo/2" => Some(Self::Indo2),
            "indo3" | "indo/3" => Some(Self::Indo3),
            "zindo1" | "zindo/1" => Some(Self::Zindo1),
            "zindo2" | "zindo/2" => Some(Self::Zindo2),
            "zindo3" | "zindo/3" => Some(Self::Zindo3),
            "indos" | "indo/s" | "zindo" | "zindos" | "zindo/s" => Some(Self::ZindoS),
            "mindo1" | "mindo/1" => Some(Self::Mindo1),
            "mindo2" | "mindo/2" => Some(Self::Mindo2),
            "mindo3" | "mindo/3" | "mindo" => Some(Self::Mindo3),
            "sindo" | "sindo1" | "sindo/1" => Some(Self::Sindo),
            "msindo" => Some(Self::Msindo),
            "mndo" => Some(Self::Mndo),
            "mndod" | "mndo/d" => Some(Self::MndoD),
            _ => None,
        }
    }
}

/// Return native methods executable by a named public API in this release.
///
/// API names are the same strings returned by [`Method::api_names`], e.g.
/// `single_point`, `gradient`, `hessian`, `excited_states`,
/// `excited_properties`, and `uv_vis_spectrum`. Call
/// [`Method::accepted_strings`] on each returned
/// method to obtain all user-facing spellings.
pub fn methods_for_api(api: &str) -> Vec<Method> {
    Method::ALL
        .into_iter()
        .filter(|method| method.api_names().contains(&api))
        .collect()
}

impl fmt::Display for Method {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for Method {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s).ok_or_else(|| format!("unknown xndo-rs method: {s}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn requested_names_are_parseable() {
        for name in [
            "cndo1", "cndo2", "cndo3", "indo", "indo1", "indo2", "indo3", "zindo1", "zindo2",
            "zindo3", "mindo1", "mindo2", "mindo3", "sindo", "msindo", "mndo", "mndod", "zindo/s",
        ] {
            assert!(Method::parse(name).is_some(), "failed to parse {name}");
        }
    }

    #[test]
    fn cndo2_and_bare_indo_are_native_single_point_methods() {
        for method in [Method::Cndo2, Method::Indo] {
            assert_eq!(method.status(), MethodStatus::Native);
            assert!(method.supports_energy());
            assert!(method.supports_uhf());
            assert!(method.api_names().contains(&"single_point"));
        }
        assert_eq!(Method::parse("indo"), Some(Method::Indo));
        assert_eq!(Method::parse("indo/s"), Some(Method::ZindoS));
    }

    #[test]
    fn noncanonical_threes_are_not_executable() {
        assert_eq!(Method::Cndo3.status(), MethodStatus::NonCanonical);
        assert_eq!(Method::Indo3.status(), MethodStatus::NonCanonical);
        assert_eq!(Method::Zindo3.status(), MethodStatus::NonCanonical);
    }
}
