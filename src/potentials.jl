# =============================================================================
# potentials.jl -- analytic Fourier components of the moire potentials
#
# Every model here is defined in real space in the paper.  This file turns
# those closed-form expressions into explicit lists of
#
#       (momentum transfer Delta_Q , Fourier amplitude)
#
# which kmodels.jl then assembles into h_QQ'(k) on the plane-wave basis.
#
# --- Real-space Fourier conventions, as printed -----------------------------
#
#   +-K valley  (B16):   V_{eta,l}(r) = sum_{G_M} e^{-eta i G_M.r} V_{eta,l,G_M}
#                        t_eta(r)     = sum_{G_M} e^{-eta i (q_1+G_M).r} t_{eta,q_1+G_M}
#
#   Gamma valley (B33):  V_{Gamma,l}(r) = sum_G V_{Gamma,l,G} e^{+i G.r}
#                        t_Gamma(r)     = sum_G t_{Gamma,G} e^{+i G.r}
#
# The two differ in the sign of the exponent; both are reproduced verbatim.
#
# --- Matrix-element rule ----------------------------------------------------
#
# (B22)/(C16) expand the *field operator*, not the state:
#
#       c†_{eta,l,r} = (1/sqrt V) sum_k sum_Q exp(-i(k-Q).r) c†_{eta,l,k-Q} .
#
# Both sides are operators, so the state created by c†_{eta,l,k-Q} carries the
# *conjugate* kernel, exp(+i(k-Q).r).  Hence, for O(r) = sum_G O_G exp(-i eta G.r),
#
#       <k,Q|O|k,Q'> = int dr e^{-i(k-Q).r} O(r) e^{+i(k-Q').r}
#                    = sum_G O_G delta(Q'-Q+eta G) = O_{eta, eta(Q-Q')} ,
#
# which is exactly the rule used by the printed momentum-space equations
# (B26)/(B27) and (C8)/(C18):
#
#       h_{Q,Q'} = O_{eta, eta(Q-Q')}                              (:paper)
#
# Reading the kernel of (B22) as the state wavefunction instead (i.e. using
# |k,Q> ∝ exp(-i(k-Q).r)) would give the opposite index ordering,
#
#       <k,Q|O|k,Q'> = O_{eta, eta(Q'-Q)}                          (:transposed)
#
# which is the plain transpose of the `:paper` matrix, i.e. the same model with
# psi -> -psi.  Since H is Hermitian, transpose = complex conjugate, so that
# variant has *identical* eigenvalues; the two are therefore physically
# indistinguishable and the comparison between them is a useful guard against a
# sign slip in the amplitude assignment.  Both are implemented and the identity
# H_paper == transpose(H_transposed) is asserted in scripts/03_selftest.jl.
# `:paper` is the default: it is the rule that (B22) actually implies.
# =============================================================================

"""one Fourier component: momentum transfer + amplitude"""
struct FourTerm
    dQ::Vector{Float64}
    c::ComplexF64
end

FourTerm(dQ, c::Real) = FourTerm(collect(float.(dQ)), ComplexF64(c))
FourTerm(dQ, c::Complex) = FourTerm(collect(float.(dQ)), ComplexF64(c))

# -----------------------------------------------------------------------------
# AA stacking, +-K valley:  (B21) + (B20), i.e. Eqs. (B17) and (B20) combined
# -----------------------------------------------------------------------------
"""
    aa_k_intralayer(geo, p, sigma; convention = :paper)

Intralayer potential (B21) for a state whose layer sign is
`sigma = (-1)^l = eta * f`.  Returns the list of Fourier components.

(B21) puts the amplitude `V exp(-sigma i psi)` on `exp(+i g_i.r)`, i.e. in the
(B16) convention `V = sum_G exp(-eta i G.r) V_G` the Fourier coefficient at
`+g_i` is `V exp(-sigma i psi)`.  Since `sigma = eta (-)^Q = eta f`, that is
exactly the `exp(-eta (-)^Q i psi)` printed in (B27).  The matrix-element rule
of the header, `h_{Q,Q'} = V_{eta,l, eta(Q-Q')}`, then assigns that same
amplitude to the entry at `Q' = Q + g_i` -- the `:paper` column.  The
`:transposed` column is the opposite index ordering, i.e. `psi -> -psi`.
"""
function aa_k_intralayer(geo::MoireGeometry, p::KParams, sigma::Int;
                         convention::Symbol = :paper)
    sgn = convention === :paper ? -1.0 : +1.0
    out = FourTerm[]
    for i in 1:3
        # +g_i  (first harmonic)
        push!(out, FourTerm(geo.g[i], p.V * exp(im * sgn * sigma * p.psi)))
        # -g_i  (first harmonic)
        push!(out, FourTerm(-geo.g[i], p.V * exp(-im * sgn * sigma * p.psi)))
        # +-g_{2i} (second harmonic, (B20): 2 V_2 sum cos = V_2 (e^{+i} + e^{-i}))
        push!(out, FourTerm(geo.g2[i], p.V2))
        push!(out, FourTerm(-geo.g2[i], p.V2))
    end
    return out
end

"""
    aa_k_interlayer(geo, p, eta)

Interlayer tunnelling (B21):  t_eta(r) = w sum_i exp(-eta i q_i.r)
                                         + w_2 sum_i exp(-eta i q_{2i}.r).
"""
function aa_k_interlayer(geo::MoireGeometry, p::KParams, eta::Int)
    out = FourTerm[]
    for i in 1:3
        push!(out, FourTerm(geo.q[i], p.w))
        push!(out, FourTerm(-geo.q[i], p.w))
        push!(out, FourTerm(geo.q2[i], p.w2))
        push!(out, FourTerm(-geo.q2[i], p.w2))
    end
    return out
end

# -----------------------------------------------------------------------------
# AB stacking, +-K valley:  Eq. (5) / (C4), momentum form (C8)
# -----------------------------------------------------------------------------
"""
    ab_k_intralayer(geo, p; convention = :paper)

Eq. (C4):  V_{eta,l}(r) = V e^{-i psi} sum_i e^{+i g_i.r} + c.c.
The potential is *layer independent* (the AB stacking has no `(-1)^l`) and has
no second harmonic (Table II fixes no V_2 for AB).
"""
function ab_k_intralayer(geo::MoireGeometry, p::KParams; convention::Symbol = :paper)
    sgn = convention === :paper ? -1.0 : +1.0
    out = FourTerm[]
    for i in 1:3
        push!(out, FourTerm(geo.g[i], p.V * exp(im * sgn * p.psi)))
        push!(out, FourTerm(-geo.g[i], p.V * exp(-im * sgn * p.psi)))
    end
    return out
end

"""
    ab_k_interlayer(geo, p)

Eq. (5)/(C4):  t_eta(r) = w sum_i exp(i eta (i-1) 2pi/3) exp(-eta i q_i.r).

In the (B16) form `sum_G exp(-eta i G.r) t_{eta,G}` this gives
`t_{eta,q_i} = w exp(i eta (i-1) 2pi/3)`.  Applying the matrix-element rule
`h_{Q,Q'} = t_{eta, eta(Q_b-Q_t)}` and using that for eta = +1 the family
`f = +1` is the bottom layer while for eta = -1 it is the top layer, the two
non-vanishing entries are

    Q' = Q - q_i   (reached from f = +1 rows) :  w exp(+i (i-1) 2pi/3)
    Q' = Q + q_i   (reached from f = -1 rows) :  w exp(-i (i-1) 2pi/3)

exactly as printed in (C8):  `w sum_i sum_{alpha=+-1} exp(i alpha (i-1) 2pi/3)
delta_{Q, Q' + alpha q_i}`.  Note the result does not depend on eta.

With AB's fitted `w = 0` (Table II) this block vanishes: the two layers
decouple, so every band is two-fold degenerate in the layer index.
"""
function ab_k_interlayer(geo::MoireGeometry, p::KParams)
    out = FourTerm[]
    for i in 1:3
        ph = exp(im * (i - 1) * 2pi / 3)          # exp(+i (i-1) 2pi/3)
        push!(out, FourTerm(-geo.q[i], p.w * ph))
        push!(out, FourTerm(geo.q[i], p.w * conj(ph)))
    end
    return out
end

# -----------------------------------------------------------------------------
# Gamma valley, AA stacking: (B38)
# -----------------------------------------------------------------------------
"""
    aa_gamma_intralayer(geo, p, gl; convention = :paper)

(B38): `V_{Gamma,l}(r) = sum_i V_Gamma exp(i (-)^l psi_Gamma) exp(i g_i.r)
                       + sum_i V_Gamma exp(-i (-)^l psi_Gamma) exp(-i g_i.r)`

with the (B33) Fourier convention `V(r) = sum_G V_G exp(+i G.r)`, so the printed
amplitudes are `V_{Gamma,+g_i} = V_Gamma exp(+i (-)^l psi_Gamma)`.
`gl = (-1)^l` is passed in (see `layer_sign_gamma`).
"""
function aa_gamma_intralayer(geo::MoireGeometry, p::GammaParams, gl::Int;
                             convention::Symbol = :paper)
    sgn = convention === :paper ? +1.0 : -1.0
    out = FourTerm[]
    for i in 1:3
        push!(out, FourTerm(geo.g[i], p.V * exp(im * sgn * gl * p.psi)))
        push!(out, FourTerm(-geo.g[i], p.V * exp(-im * sgn * gl * p.psi)))
    end
    return out
end

# -----------------------------------------------------------------------------
# Gamma valley, AB stacking: (C13)
# -----------------------------------------------------------------------------
"""
    ab_gamma_intralayer(geo, p; convention = :paper)

(C13):  V_{Gamma,l}(r) = [ sum_i V_Gamma e^{+i psi_Gamma} e^{+i g_i.r}
                         + sum_i V_Gamma e^{-i psi_Gamma} e^{-i g_i.r} ] s_0

Layer *independent*, and diagonal in spin (hence every band is at least two-fold
degenerate; see the text after (C18)).
"""
function ab_gamma_intralayer(geo::MoireGeometry, p::GammaParams; convention::Symbol = :paper)
    sgn = convention === :paper ? +1.0 : -1.0
    out = FourTerm[]
    for i in 1:3
        push!(out, FourTerm(geo.g[i], p.V * exp(im * sgn * p.psi)))
        push!(out, FourTerm(-geo.g[i], p.V * exp(-im * sgn * p.psi)))
    end
    return out
end

"""
    gamma_interlayer(p)

`t_Gamma(r) = w_Gamma s_0` ((B38) for AA, (C11) for AB): a pure `G = 0`
transfer, diagonal in spin.  Returns the single component.
"""
gamma_interlayer(p::GammaParams) = [FourTerm([0.0, 0.0], p.w)]

# The public API passes a ParamSet; these thin methods unwrap it so that the
# constructors above only ever see a bare payload.
aa_k_intralayer(geo::MoireGeometry, ps::ParamSet, sigma::Int; kw...) =
    aa_k_intralayer(geo, payload(ps), sigma; kw...)
aa_k_interlayer(geo::MoireGeometry, ps::ParamSet, eta::Int) =
    aa_k_interlayer(geo, payload(ps), eta)
ab_k_intralayer(geo::MoireGeometry, ps::ParamSet; kw...) =
    ab_k_intralayer(geo, payload(ps); kw...)
ab_k_interlayer(geo::MoireGeometry, ps::ParamSet) = ab_k_interlayer(geo, payload(ps))
aa_gamma_intralayer(geo::MoireGeometry, ps::ParamSet, gl::Int; kw...) =
    aa_gamma_intralayer(geo, payload(ps), gl; kw...)
ab_gamma_intralayer(geo::MoireGeometry, ps::ParamSet; kw...) =
    ab_gamma_intralayer(geo, payload(ps); kw...)
gamma_interlayer(ps::ParamSet) = gamma_interlayer(payload(ps))

# -----------------------------------------------------------------------------
# integer key increments
# -----------------------------------------------------------------------------
"""
    dkey(geo, v, scale) -> (dn1, dn2)

Integer key increment for a momentum transfer `v`, for a basis whose states
satisfy `scale * Q = n1 b_{M,1} + n2 b_{M,2}` (scale = 3 for the +-K valley
basis of (B23), scale = 1 for the Gamma-valley basis).  Asserts that `v` really
lands on the key grid -- i.e. that the Fourier component list and the basis are
mutually consistent.
"""
function dkey(geo::MoireGeometry, v, scale::Int; tol::Float64 = 1e-7)
    c1, c2 = lattice_coeffs(geo, v)
    n1, n2 = c1 * scale, c2 * scale
    @assert abs(n1 - round(n1)) < tol && abs(n2 - round(n2)) < tol "momentum transfer $(v) is off the key grid"
    return (round(Int, n1), round(Int, n2))
end

"""convert a Fourier-component list into (integer key increment, amplitude) pairs"""
function to_dkeys(geo::MoireGeometry, terms::Vector{FourTerm}, scale::Int)
    out = Tuple{NTuple{2,Int},ComplexF64}[]
    for t in terms
        push!(out, (dkey(geo, t.dQ, scale), t.c))
    end
    # combine duplicates (e.g. the second-harmonic list can degenerate)
    acc = Dict{NTuple{2,Int},ComplexF64}()
    for (k, c) in out
        acc[k] = get(acc, k, zero(ComplexF64)) + c
    end
    # produce real tuples, not the Pairs that `collect` on a Dict would give
    return Tuple{NTuple{2,Int},ComplexF64}[(k, c) for (k, c) in acc]
end
