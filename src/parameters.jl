# =============================================================================
#   model parameters from papers
#
# All energies in meV, masses in units of the free electron mass m_e,
# momenta in Angstrom^-1.
#
# A parameter set is a `ParamSet`: a physics payload (`KParams` / `GammaParams`)
# wrapped together with a `ParamRef` that records WHERE the numbers came from and
# WHICH geometry they are valid at.  The reason is that Table I and Table II are
# fits to the DFT bands **at theta = 3.89 deg, for MoTe2 only** -- they are not
# material constants.  The FH+SH set in particular uses its second harmonics to
# encode lattice relaxation, which grows with the moire period: the DFT top
# valence bandwidth is 12.8 meV at 3.89 deg and 133.3 meV at 9.43 deg (Table V).
# Pairing a parameter set with a geometry it was not fitted at is a real physics
# error that does not announce itself, so `check_parameters` (called at the top of
# `hamiltonian` and `transcribe_hamiltonian`) turns it into a loud one.
#
# Note on `const`: these bindings stay `const`.  For an immutable record that is
# idiomatic and does NOT restrict generality -- any number of further `ParamSet`s
# (other materials, angles, harmonics) can be defined alongside these; `const`
# only forbids rebinding the same name.  Dropping it would be worse: a non-const
# global is typed `Any` at each use site.
# =============================================================================

# -----------------------------------------------------------------------------
# provenance and validity
# -----------------------------------------------------------------------------
"""how a parameter set was obtained"""
@enum Provenance FITTED CALIBRATED PLACEHOLDER EXPLORATORY

"""
    ParamRef(; material, a0, theta_deg, stacking, valley, harmonics,
               provenance, source, calibrated_at)

Where a parameter set came from, and the geometry it is valid at.

`a0` and `theta_deg` use `NaN` for "this record makes no claim about its
geometry".  That is what `ParamRef()` -- and hence any hand-built bare
`KParams` / `GammaParams` -- produces, so a record that says nothing is exempt
from `check_parameters` and hand-built parameter sets stay unconstrained by
default.  `provenance = EXPLORATORY` is the explicit form of the same statement.
"""
struct ParamRef
    material::Symbol       # :MoTe2, :unknown, ...
    a0::Float64            # valid at this monolayer lattice constant (Ang); NaN = no claim
    theta_deg::Float64     # valid at this twist angle (deg); NaN = no claim
    stacking::Symbol       # :AA, :AB, :unknown
    valley::Symbol         # :K, :Gamma, :unknown
    harmonics::Symbol      # :FH, :FHSH, :none
    provenance::Provenance
    source::String         # which paper / table / rule
    calibrated_at::Union{Nothing,Tuple{Float64,Float64}}   # (a0, theta_deg) for CALIBRATED
end

ParamRef(; material::Symbol = :unknown, a0::Real = NaN, theta_deg::Real = NaN,
         stacking::Symbol = :unknown, valley::Symbol = :unknown,
         harmonics::Symbol = :none, provenance::Provenance = EXPLORATORY,
         source::AbstractString = "", calibrated_at = nothing) =
    ParamRef(material, float(a0), float(theta_deg), stacking, valley, harmonics,
             provenance, String(source), calibrated_at)

# -----------------------------------------------------------------------------
# materials -- `a0` is geometry, not a fitted number, so it gets its own record
# -----------------------------------------------------------------------------
"""
    Material(name, a0, display)

A monolayer lattice constant.  Nothing in these continuum models is
material-specific except the numbers fitted to it, and those live in a
`ParamSet`; a `Material` exists so that `a0` stops being a bare literal.

To add a material, add one of these next to its own `ParamSet`s -- there is no
registry to update.  Note that PRB 109, 205121 (2024) fits **only** MoTe2, so a
new material needs its own V, psi, w, ... from whatever source you are using;
its geometry needs nothing but `a0`.
"""
struct Material
    name::Symbol
    a0::Float64
    display::String
end

"""monolayer lattice constant of MoTe2.  Experimental a = 3.519 Ang (App. A1);
DFT-D2 gives 3.518.  The paper's fits use a_0 = 3.52 Ang."""
const A0_MOTE2 = A0_DEFAULT

const MOTE2 = Material(:MoTe2, A0_MOTE2, "MoTe2")

"""`moire_geometry(MOTE2)` == `moire_geometry(a0 = 3.52, ...)`; same for any Material."""
moire_geometry(mat::Material; theta_deg::Real = THETA_DEFAULT_DEG, Nshell::Int = 5) =
    moire_geometry(; a0 = mat.a0, theta_deg = theta_deg, Nshell = Nshell)

# -----------------------------------------------------------------------------
# the wrapper: metadata + payload
# -----------------------------------------------------------------------------
"""supertype of the two physics payloads"""
abstract type AbstractParams end

"""
    ParamSet{P<:AbstractParams}

A parameter set: `params` is the physics payload, `ref` says where it came from
and which geometry it is valid at.

The public API takes a `ParamSet`; the assembly kernels (`ham_AA_K`, the
potentials in potentials.jl, the transcriptions in transcribe.jl) still take a
bare payload and dispatch on its concrete type, so `payload` is unwrapped at the
boundary -- see `hamiltonian` in kmodels.jl.
"""
struct ParamSet{P<:AbstractParams}
    ref::ParamRef
    params::P
end

"""the physics payload, wrapped or not"""
payload(ps::ParamSet) = ps.params
payload(p::AbstractParams) = p

"""the provenance record; a bare payload makes no claim"""
paramref(ps::ParamSet) = ps.ref
paramref(::AbstractParams) = ParamRef()

label(ps::ParamSet) = ps.params.label
label(p::AbstractParams) = p.label

"""
    exploratory(ps::ParamSet) -> ParamSet

The same payload carrying a provenance record that makes no claim about its
geometry, so `check_parameters` will not fire.  Use this when you *deliberately*
want to evaluate paper parameters at another material or angle.
"""
exploratory(ps::ParamSet) = ParamSet(
    ParamRef(material = paramref(ps).material, stacking = paramref(ps).stacking,
             valley = paramref(ps).valley, harmonics = paramref(ps).harmonics,
             provenance = EXPLORATORY, source = paramref(ps).source), ps.params)

# -----------------------------------------------------------------------------
# +/-K valley: (B15)/(B17)/(B20)/(B21) for AA and (C4)/(C8) for AB
# -----------------------------------------------------------------------------
struct KParams <: AbstractParams
    mstar::Float64     # effective mass m* (m_e), same in both layers and valleys
    V::Float64         # first-harmonic potential amplitude
    psi::Float64       # potential phase (radians)
    w::Float64         # first-harmonic interlayer tunnelling
    V2::Float64        # second harmonic potential amplitude   ((B20); AA only)
    w2::Float64        # second harmonic tunnelling amplitude  ((B20); AA only)
    label::String
end

KParams(; mstar, V, psi, w, V2 = 0.0, w2 = 0.0, label = "") =
    KParams(float(mstar), float(V), float(psi), float(w), float(V2), float(w2), label)

# The five sets below are records of PRB 109, 205121 (2024), every one of which is
# a single-angle fit.  Their `a0` / `theta_deg` name the geometry through the
# package defaults (`A0_MOTE2` == `A0_DEFAULT`, `THETA_DEFAULT_DEG`), which were
# taken from that paper -- so the two coincide today, but the records are pinned
# to the defaults rather than to the paper.  Changing either default therefore
# moves the geometry these five claim, and Part G of the self-test asserts both
# defaults are still 3.52 / 3.89 so that it cannot happen unnoticed.

"""Table I -- AA stacking, +-K valley, first harmonics only (V2 = w2 = 0)"""
const AA_K_FH = ParamSet(
    ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG,
             stacking = :AA, valley = :K, harmonics = :FH, provenance = FITTED,
             source = "PRB 109, 205121 (2024), Table I; fit to the DFT bands at theta = 3.89 deg"),
    KParams(mstar = 0.60, V = 16.5, psi = deg2rad(-105.9), w = -18.8,
            V2 = 0.0, w2 = 0.0, label = "AA +/-K FH"))

"""Table I -- AA stacking, +-K valley, first + second harmonics"""
const AA_K_FHSH = ParamSet(
    ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG,
             stacking = :AA, valley = :K, harmonics = :FHSH, provenance = FITTED,
             source = "PRB 109, 205121 (2024), Table I; fit to the DFT bands at theta = 3.89 deg.  The " *
                      "second harmonics stand in for lattice relaxation, which is " *
                      "what makes this set the most angle-sensitive of the five"),
    KParams(mstar = 0.62, V = 7.94, psi = deg2rad(-88.43), w = -10.77,
            V2 = 20.00, w2 = 10.21, label = "AA +/-K FH+SH"))

"""Table II -- AB stacking, +-K valley.  The AB +/-K model (C8) has no SH terms."""
const AB_K = ParamSet(
    ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG,
             stacking = :AB, valley = :K, harmonics = :FH, provenance = FITTED,
             source = "PRB 109, 205121 (2024), Table II; fit to the DFT bands at " *
                      "theta = 3.89 deg"),
    KParams(mstar = 0.62, V = 53.0, psi = deg2rad(-56.0), w = 0.0,
            V2 = 0.0, w2 = 0.0, label = "AB +/-K"))

# -----------------------------------------------------------------------------
# Gamma valley: (B32)-(B38) for AA and (C11)-(C18) for AB
# -----------------------------------------------------------------------------
struct GammaParams <: AbstractParams
    mstar::Float64     # m_Gamma* (m_e)
    V::Float64         # V_Gamma
    psi::Float64       # psi_Gamma (radians)
    w::Float64         # w_Gamma
    E0::Float64        # E_Gamma: free bias, NOT tabulated (see note below)
    label::String
end

GammaParams(; mstar, V, psi, w, E0 = 0.0, label = "") =
    GammaParams(float(mstar), float(V), float(psi), float(w), float(E0), label)

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
const AB_GAMMA = ParamSet(
    ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG,
             stacking = :AB, valley = :Gamma, harmonics = :none,
             provenance = CALIBRATED,
             source = "PRB 109, 205121 (2024), Table II for m_Gamma*, V_Gamma and " *
                      "w_Gamma; E_Gamma is not tabulated anywhere and is calibrated " *
                      "here by the Fig. 4(c) '30 meV below the +-K VBM' rule",
             calibrated_at = (A0_MOTE2, THETA_DEFAULT_DEG)),
    GammaParams(mstar = 10.0, V = 72.0, psi = 0.0, w = 300.0, E0 = -640.1296,
                label = "AB Gamma"))

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
const AA_GAMMA = ParamSet(
    ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG,
             stacking = :AA, valley = :Gamma, harmonics = :none,
             provenance = PLACEHOLDER,
             source = "no fit exists in PRB 109, 205121 (2024) -- App. B 3 fits the " *
                      "+-K model alone.  The AB values are reused (at psi_Gamma = 0 " *
                      "the AA and AB potentials coincide) and E_Gamma is a placeholder"),
    GammaParams(mstar = 10.0, V = 72.0, psi = 0.0, w = 300.0, E0 = -817.3836,
                label = "AA Gamma (unfitted)"))

# -----------------------------------------------------------------------------
# the paper's fits, indexed by (stacking, valley, harmonics)
# -----------------------------------------------------------------------------
"""
    PAPER_FITS

Lookup index of the recorded fits, derived from the constants above -- a table,
not a second source of truth.  Keys are `(stacking, valley, harmonics)`.
"""
const PAPER_FITS = Dict{Tuple{Symbol,Symbol,Symbol},ParamSet}(
    (:AA, :K,     :FH)   => AA_K_FH,
    (:AA, :K,     :FHSH) => AA_K_FHSH,
    (:AB, :K,     :FH)   => AB_K,
    (:AA, :Gamma, :none) => AA_GAMMA,
    (:AB, :Gamma, :none) => AB_GAMMA)

"""
    paper_fit(stacking, valley, harmonics = :none) -> ParamSet

Recorded fit, e.g. `paper_fit(:AB, :Gamma)`.  Throws, listing the available keys,
if there is no such fit.
"""
function paper_fit(stacking::Symbol, valley::Symbol, harmonics::Symbol = :none)
    key = (stacking, valley, harmonics)
    haskey(PAPER_FITS, key) ||
        error("no recorded fit for $key; available: " *
              join(string.(sort(collect(keys(PAPER_FITS)))), ", "))
    return PAPER_FITS[key]
end

# -----------------------------------------------------------------------------
# validity guard: a parameter set meets a geometry only where it was fitted
# -----------------------------------------------------------------------------
const ENFORCE_PARAM_CHECK = Ref(true)

"""
    enforce_param_checks!(flag::Bool)

Session-wide switch for `check_parameters` (on by default).  Prefer the
per-record `exploratory(ps)`, which says the same thing locally.
"""
enforce_param_checks!(flag::Bool) = (ENFORCE_PARAM_CHECK[] = flag)

"""which (stacking, valley) each model symbol stands for"""
const MODEL_REF = Dict(:AA_K => (:AA, :K), :AB_K => (:AB, :K),
                       :AA_Gamma => (:AA, :Gamma), :AB_Gamma => (:AB, :Gamma))

"""
    check_parameters(geo, p; model = nothing, atol_a0 = 1e-4, atol_theta_deg = 1e-4)

Throw if parameter set `p` is used at a geometry it was not fitted at, or with a
Hamiltonian it does not belong to.

Table I / Table II are fits to the DFT bands at theta = 3.89 deg for MoTe2 -- the
DFT top valence bandwidth moves from 12.8 meV at 3.89 deg to 133.3 meV at 9.43 deg
(Table V) -- so a mismatch is a physics error, not a rounding question.  It is
silent otherwise: the bands still come out smooth and plausible.

A record that makes no claim (`a0` or `theta_deg` NaN, or provenance
`EXPLORATORY`) is skipped, so hand-built parameter sets are unaffected.  To
deliberately cross a material or angle, wrap with `exploratory(ps)`, or switch
the guard off for the session with `enforce_param_checks!(false)`.
"""
function check_parameters(geo::MoireGeometry, p; model::Union{Nothing,Symbol} = nothing,
                          atol_a0::Float64 = 1e-4, atol_theta_deg::Float64 = 1e-4)
    ENFORCE_PARAM_CHECK[] || return true
    ref = paramref(p)
    (isfinite(ref.a0) && isfinite(ref.theta_deg) && ref.provenance !== EXPLORATORY) ||
        return true

    theta_deg = rad2deg(geo.theta)
    if !isapprox(geo.a0, ref.a0; atol = atol_a0)
        error("parameter set '$(label(p))' is valid at a0 = $(ref.a0) Ang " *
              "($(ref.material)), but the geometry has a0 = $(geo.a0) Ang.  " *
              "Build the geometry from the matching Material, or wrap the parameters " *
              "in exploratory(ps) if the mismatch is intended.")
    end
    if !isapprox(theta_deg, ref.theta_deg; atol = atol_theta_deg)
        error("parameter set '$(label(p))' is valid at theta = $(ref.theta_deg) deg, " *
              "but the geometry is at theta = $(theta_deg) deg.  The tabulated fits " *
              "are all made at a single angle; refit for this angle, or wrap the " *
              "parameters in exploratory(ps) if the mismatch is intended.")
    end
    if model !== nothing && haskey(MODEL_REF, model)
        st, va = MODEL_REF[model]
        (ref.stacking === st && ref.valley === va) ||
            error("parameter set '$(label(p))' is for $(ref.stacking)/$(ref.valley) " *
                  "but model = $(model) is $(st)/$(va).")
    end
    return true
end

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

# the public API passes a ParamSet; unwrap to the payload for the formulas above
hbar_omega(geo::MoireGeometry, ps::ParamSet) = hbar_omega(geo, payload(ps))
harmonic_stiffness(geo::MoireGeometry, ps::ParamSet) = harmonic_stiffness(geo, payload(ps))
anharmonic_shift(geo::MoireGeometry, ps::ParamSet) = anharmonic_shift(geo, payload(ps))

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
