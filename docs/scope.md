# Scope

Version 0.2.4 provides molecular, non-periodic calculations for the native
methods listed in `README.md`. Every native ground-state engine exposes an
analytic gradient and Hessian with RHF/UHF references. The ZINDO/S s/p branch
adds spin-adapted singlet/triplet RHF-CIS and spin-orbital UCIS properties,
excitation derivatives, and absolute-state derivatives.

The following are outside this release:

- PM3 and all PM3 parameter data;
- D3, H4, HX, and other empirical corrections;
- Sparkles, point charges, and capped-bond pseudoatoms;
- periodic boundaries and solvation;
- transition-metal 9-AO ZINDO/S;
- transition-metal 9-AO ZINDO/S state derivatives;
- excited-state geometry optimization and nonadiabatic derivative couplings;
- registered methods without a complete redistributable parameter set.

Unsupported methods fail explicitly. Values from a related method are never
silently substituted.

Open-shell diatomics, including OH, retain UHF energies, analytic gradients,
and vertical UCIS spectra. Their state-specific UHF Hessians and UCIS
derivatives are not unique on a degenerate electronic surface unless a
state-averaged or diabatic model is selected; xndo-rs therefore rejects those
derivative requests explicitly.
Requested CIS/UCIS roots closer than 0.0001 eV are treated the same way:
vertical properties are retained, but a branch-dependent state derivative is
not returned as if it were unique.
