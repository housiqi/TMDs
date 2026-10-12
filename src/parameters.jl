# =============================================================================
# parameters.jl -- model parameters, transcribed from the paper's tables
#
#   Table I  : AA-stacked tMoTe2, +/-K valley, Eq. (3) / (B27)
#   Table II : AB-stacked tMoTe2, +/-K valley Eq. (5)/(C8) and Gamma valley
#              Eq. (7)/(C18)
#
# All energies in meV, masses in units of the free electron mass m_e,
# momenta in Angstrom^-1.
# =============================================================================

"""monolayer lattice constant of MoTe2.  Experimental a = 3.519 Ang (App. A1);
DFT-D2 gives 3.518.  The paper's fits use a_0 = 3.52 Ang."""
const A0_MOTE2 = 3.52

const THETA_DEG = 3.89        # the twist angle used for every fit in the paper

# -----------------------------------------------------------------------------
# +/-K valley: (B15)/(B17)/(B20)/(B21) for AA and (C4)/(C8) for AB
# -----------------------------------------------------------------------------
struct KParams
    mstar::Float64     # effective mass m* (m_e), same in both layers and valleys
    V::Float64         # first-harmonic potential amplitude
    psi::Float64       # potential phase (radians)
    w::Float64         # first-harmonic interlayer tunnelling
    V2::Float64        # second harmonic potential amplitude   ((B20); AA only)
    w2::Float64        # second harmonic tunnelling amplitude  ((B20); AA only)
    label::String
end

"""Table I -- AA stacking, +-K valley, first harmonics only (V2 = w2 = 0)"""
const AA_K_FH = KParams(0.60, 16.5, deg2rad(-105.9), -18.8, 0.0, 0.0, "AA +/-K FH")

"""Table I -- AA stacking, +-K valley, first + second harmonics"""
const AA_K_FHSH = KParams(0.62, 7.94, deg2rad(-88.43), -10.77, 20.00, 10.21, "AA +/-K FH+SH")

"""Table II -- AB stacking, +-K valley.  The AB +/-K model (C8) has no SH terms."""
const AB_K = KParams(0.62, 53.0, deg2rad(-56.0), 0.0, 0.0, 0.0, "AB +/-K")

# -----------------------------------------------------------------------------
# Gamma valley: (B32)-(B38) for AA and (C11)-(C18) for AB
# -----------------------------------------------------------------------------
struct GammaParams
    mstar::Float64     # m_Gamma* (m_e)
    V::Float64         # V_Gamma
    psi::Float64       # psi_Gamma (radians)
    w::Float64         # w_Gamma
    E0::Float64        # E_Gamma: free bias, NOT tabulated (see note below)
    label::String
end

"""
Table II -- AB stacking, Gamma valley.

`E0 = E_Gamma`.  **It is not given in any table of the paper**, and cannot be
derived from the tables: (C19) is `h = hbar^2 grad^2/2m_Gamma* + w_Gamma +
E_Gamma + V_Gamma(r)`, in which `E_Gamma` is a rigid shift of the whole
Gamma-valley manifold.  App. B 3 says as much -- the zeroth harmonic "must be
included to compare the energies of the in-active Gamma valley and active +-K
valleys" -- but never tabulates it, because each valley is fitted separately.
It is therefore a free parameter, and the default below is *calibrated* (not
fitted) by the paper's own statement (Sec. III.C; Fig. 4(c) caption) that the
first pair of ultra-flat Gamma bands lies 30 meV below the +-K valence band
maximum:

    E0 = (VBM of the +-K model) - 30 meV - (top Gamma band at E0 = 0)
       = 66.053 - 30 - 676.183 = -640.130 meV        (theta = 3.89 deg)

Both VBM and flat-band energies are the maxima along the high-symmetry path at
this angle; `E0_relative_to_vbm` in bands.jl recomputes this for other angles,
and scripts/03_selftest.jl asserts the resulting position.  At psi_Gamma = 0
only the *offset* depends on E0, never the band structure or the 54.9 meV
splitting of the two flat bands.
"""
const AB_GAMMA = GammaParams(10.0, 72.0, 0.0, 300.0, -640.1296, "AB Gamma")

"""
AA stacking, Gamma valley.

**No fitted parameters exist for this model.**  App. B 3 states explicitly that
the fits used the +-K valley model alone (Eq. (B27)) because the Gamma valley
lies below the 6th valence band for theta <= 5.09 deg, so Table I contains no
Gamma parameters.

We therefore expose all of them and use the AB values as the default: at
`psi_Gamma = 0` (Table II) the AA and AB Gamma-valley potentials (B38) and
(C13) are the *same* expression, so this is the least arbitrary choice
available.  The offset is separate and is set below by the paper's own
qualitative statement for this stacking -- the Gamma valley is out of the
plotted window, below the 6th valence band -- i.e. 160 meV under the +-K
valence band maximum.  That number is a placeholder, not a fit.
"""
const AA_GAMMA = GammaParams(10.0, 72.0, 0.0, 300.0, -817.3836, "AA Gamma (unfitted)")

# -----------------------------------------------------------------------------
# derived scales
# -----------------------------------------------------------------------------

"""hbar^2/(2 m*) in meV*Angstrom^2 -- the coefficient of |k-Q|^2."""
kinetic_prefactor(mstar::Real) = HBAR2_2ME / mstar

"""
    hbar_omega(geo, gamma) -> meV

Eq. (C25):  hbar*omega = hbar*sqrt(3 V_Gamma |b_{M,1}|^2 / m_Gamma*).

Squaring, `hbar^2 * 3 V_Gamma |b|^2 / m_Gamma*` = `6 V_Gamma * kinetic_prefactor`,
which with the Table II values at theta = 3.89 deg evaluates to 56.76 meV --
the number quoted in (C25).  Asserted in scripts/03_selftest.jl.
"""
function hbar_omega(geo::MoireGeometry, p::GammaParams)
    return sqrt(6 * p.V * kinetic_prefactor(p.mstar) * norm(geo.bM1)^2)
end

"""Eq. (C22)/(C23): 2 V_Gamma sum_i cos(g_i.r) = 6V_Gamma - (3/2)V_Gamma|b_{M,1}|^2 r^2 + O(r^4)"""
function harmonic_stiffness(geo::MoireGeometry, p::GammaParams)
    return 1.5 * p.V * norm(geo.bM1)^2
end

"""
    anharmonic_shift(geo, p) -> meV

First-order perturbative correction to the `(0,0) -> (1,0),(0,1)` splitting of
the Gamma-valley harmonic oscillator (C23) from the *quartic* term of the exact
lattice potential.

Expanding (C20) exactly, the potential that the oscillator (C23) approximates is

    -V_Gamma(r) = -6V_Gamma + (3/2)V_Gamma|b_{M,1}|^2 r^2
                           - (3/32)V_Gamma|b_{M,1}|^4 r^4 + O(r^6),

so the perturbation is `W = -(3/32)V_Gamma|b|^4 r^4`.  With the 2D oscillator
length `l^2 = hbar/(m_Gamma* omega)` the needed moments are
`<r^4>_(0,0) = 2 l^4` and `<r^4>_avg[(1,0),(0,1)] = 6 l^4`, hence

    dE_1 - dE_0 = -(3/8) V_Gamma |b|^4 l^4 = -(1/4) * hbar^2 |b_{M,1}|^2 / (2 m_Gamma*),

a *negative* shift that is **independent of V_Gamma** (the `l^4` cancels it).
At theta = 3.89 deg it is -1.865 meV, so the exact lattice splitting should be
`hbar*omega - 1.865 = 54.90 meV` rather than the `56.76 meV` of (C25).

That is the whole content of the difference from the paper's quoted 56.74 meV:
(C24)/(C25) solve the *harmonic* problem, and App. C 3 says so in as many words
-- the flat-band energies "can be approximated by an array of decoupled harmonic
oscillators".  Both numbers are checked against the numerical diagonalisation in
scripts/03_selftest.jl.
"""
anharmonic_shift(geo::MoireGeometry, p::GammaParams) =
    -kinetic_prefactor(p.mstar) * norm(geo.bM1)^2 / 4

# -----------------------------------------------------------------------------
# layer-sign conventions -- the paper uses OPPOSITE labels in App. B and App. C
# -----------------------------------------------------------------------------
"""
    layer_sign_kvalley(f, eta) -> +-1

The layer sign (-1)^l of Eqs. (B19)/(B26), expressed through the Q family.

(B17) defines `l = t, b` respectively as `l = 0, 1` for `(-)^l`, i.e. the TOP
layer carries +1.  (B23) assigns `Q in {G_M + eta q_1}` to layer `l` with
`(-1)^l = +1`, and (B26) uses `eta (-)^Q = (-1)^{l_Q^eta}`; hence with
`f = +1 <-> Q in {G_M + q_1}` and `f = -1 <-> Q in {G_M - q_1}`:

    (-)^Q = f         and        layer sign sigma(Q,eta) = eta * f.

This is exactly the combination that appears as `(eps/2) eta (-)^Q` in (B26).
"""
layer_sign_kvalley(f::Int, eta::Int) = eta * f

"""
    layer_sign_gamma(l) -> +-1

Layer sign `(-1)^l` for the Gamma-valley appendices.

(B36) reads `V_{Gamma,b,+-g_i} = V_Gamma exp(+-i psi_Gamma)` together with
(B38) `V_{Gamma,l} = sum_i V_Gamma exp(i (-)^l psi_Gamma) exp(i g_i.r) + ...`,
which fixes `(-)^b = +1`, `(-)^t = -1` -- the OPPOSITE labelling from (B17)
used in App. B.  The two appendices therefore disagree about which layer gets
`+eps/2`; this is only a naming choice (equivalent to eps -> -eps) and is
irrelevant to the band structures at eps = 0, where every result here is
computed.  At psi_Gamma = 0 (Table II) it is irrelevant even for the potential.

We keep the printed convention of App. C for the Gamma-valley models.
"""
layer_sign_gamma(l::Symbol) = l === :b ? +1 : -1
