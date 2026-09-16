# Scope

Version 0.3.0 provides molecular, non-periodic calculations for the native
methods listed in `README.md`. Every native ground-state engine exposes an
analytic gradient and Hessian with RHF/UHF references, orbital energies with the
HOMO/LUMO pair, and Molden wavefunction output. The ZINDO/S s/p branch adds
spin-adapted singlet/triplet RHF-CIS and spin-orbital UCIS properties,
excitation derivatives, and absolute-state derivatives.

`docs/v0.3.0-delivery.md` records what this release delivers against what was
planned, and states the known problems in it.

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
