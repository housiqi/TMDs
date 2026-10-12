# =============================================================================
# kmodels.jl -- assemble h_QQ'(k) on the plane-wave basis
#
# All four models of the paper reduce to diagonalising a plane-wave matrix
#     H_QQ'(k) = <k,Q| H |k,Q'>
# in the basis (B22) / (C16):   |k,Q> ∝ exp(-i (k-Q).r).
#
# Three ingredients are summed, following (B26)/(B27) for the +-K valley and
# (C17)/(C18) for the Gamma valley:
#
#   1. diagonal       -hbar^2 (k-Q)^2 / (2 m*) + (eps/2) * layer sign
#   2. intralayer     sum over the Fourier components of V(r)   (potentials.jl)
#   3. interlayer     sum over the Fourier components of t(r)   (potentials.jl)
#
# The Fourier amplitudes come from the *real-space* expressions (B21)/(B38)/
# (C4)/(C13); the momentum-space forms printed in (B26)/(B27)/(C8)/(C18) are
# reproduced verbatim in transcribe.jl and compared against this file in the
# self-test.  See the header of potentials.jl for the index-ordering caveat.
# =============================================================================

# -----------------------------------------------------------------------------
# AA stacking, +-K valley:  (B15) + (B21) -> (B27)
# -----------------------------------------------------------------------------
"""
    ham_AA_K(geo, bas, p, k; eta = 1, eps = 0.0, convention = :paper) -> Matrix

AA-stacked +-K-valley Hamiltonian of Eq. (B27) (Eq. (3) of the main text) at
momentum `k`, on the plane-wave set of (B23).  `eta = +-1` selects the valley;
`eps` is the displacement field of (B19).

`convention = :paper` uses the index ordering that (B22) implies and (B27)
prints; `:transposed` uses the opposite ordering, which is the same model with
`psi -> -psi`.  The two matrices are plain transposes of one another and
therefore have identical spectra -- see the header of potentials.jl.
"""
function ham_AA_K(geo::MoireGeometry, bas::KValleyBasis, p::KParams, k;
                  eta::Int = 1, eps::Real = 0.0, convention::Symbol = :paper)
    n = length(bas.Qs)
    H = zeros(ComplexF64, n, n)
    kin = kinetic_prefactor(p.mstar)

    # intralayer terms depend on the layer sign sigma = eta*f, i.e. on f
    Vterms = Dict(
        +1 => to_dkeys(geo, aa_k_intralayer(geo, p, layer_sign_kvalley(+1, eta);
                                            convention = convention), 3),
        -1 => to_dkeys(geo, aa_k_intralayer(geo, p, layer_sign_kvalley(-1, eta);
                                            convention = convention), 3))
    Tterms = to_dkeys(geo, aa_k_interlayer(geo, p, eta), 3)

    for a in 1:n
        ka = bas.key[a]
        f = bas.f[a]
        sigma = layer_sign_kvalley(f, eta)
        dk = k .- bas.Qs[a]
        H[a, a] += -kin * dot(dk, dk) + eps / 2 * sigma

        for (dkey_, c) in Vterms[f]
            b = lookup(bas, key_add(ka, dkey_))
            b == 0 && continue
            bas.f[b] == f || continue          # intralayer: family is conserved
            H[a, b] += c
        end
        for (dkey_, c) in Tterms
            b = lookup(bas, key_add(ka, dkey_))
            b == 0 && continue
            H[a, b] += c
        end
    end
    return H
end

# -----------------------------------------------------------------------------
# AB stacking, +-K valley:  Eq. (5) -> (C8)
# -----------------------------------------------------------------------------
"""
    ham_AB_K(geo, bas, p, k; eta = 1, eps = 0.0, convention = :paper)

AB-stacked +-K-valley Hamiltonian of Eq. (C8) (Eq. (5) of the main text).

Differences from AA: the intralayer potential is layer *independent* (no
`(-1)^l`), the tunnelling carries the `exp(i eta (i-1) 2pi/3)` phases, and there
are no second harmonics.  With the fitted `w = 0` the layer blocks decouple.
"""
function ham_AB_K(geo::MoireGeometry, bas::KValleyBasis, p::KParams, k;
                  eta::Int = 1, eps::Real = 0.0, convention::Symbol = :paper)
    n = length(bas.Qs)
    H = zeros(ComplexF64, n, n)
    kin = kinetic_prefactor(p.mstar)

    Vterms = to_dkeys(geo, ab_k_intralayer(geo, p; convention = convention), 3)
    Tterms = to_dkeys(geo, ab_k_interlayer(geo, p), 3)

    for a in 1:n
        ka = bas.key[a]
        sigma = layer_sign_kvalley(bas.f[a], eta)
        dk = k .- bas.Qs[a]
        H[a, a] += -kin * dot(dk, dk) + eps / 2 * sigma

        for (dkey_, c) in Vterms
            b = lookup(bas, key_add(ka, dkey_))
            b == 0 && continue
            bas.f[b] == bas.f[a] || continue   # intralayer: family is conserved
            H[a, b] += c
        end
        for (dkey_, c) in Tterms
            b = lookup(bas, key_add(ka, dkey_))
            b == 0 && continue
            H[a, b] += c
        end
    end
    return H
end

# -----------------------------------------------------------------------------
# Gamma valley: basis = (G_M, layer, spin)
# -----------------------------------------------------------------------------
"""
    GammaBasis

Plane-wave basis for the Gamma-valley models, ordered as
`((G_1,b,↑),(G_1,b,↓),(G_1,t,↑),(G_1,t,↓),(G_2,b,↑),...)` -- the (G_M, layer)
index runs slow and the spin index fast (`gamma_with_spin` is `kron(H, s_0)`).
Every Gamma-valley term is proportional to `s_0`, so the spectrum is exactly
two-fold spin degenerate.  `n` is the *full* dimension, i.e. 4 per G_M.
"""
struct GammaBasis
    G::GammaValleyBasis
    n::Int          # 4 * length(G.Qs)
end

gammabasis(geo::MoireGeometry) = GammaBasis(gammavalley_basis(geo), 4 * length(geo.Gs))

@inline gidx(bas::GammaBasis, gi::Int, li::Int) = (gi - 1) * 2 + li   # li = 1 -> b, 2 -> t

"""momentum-space Hamiltonian of the Gamma-valley models, spinless 2x2-layer form"""
function ham_gamma_spinless(geo::MoireGeometry, bas::GammaBasis, p::GammaParams,
                            intraterms::Vector{Tuple{NTuple{2,Int},ComplexF64}}, k;
                            eps::Real = 0.0)
    nG = length(bas.G.Qs)
    H = zeros(ComplexF64, 2nG, 2nG)
    kin = kinetic_prefactor(p.mstar)
    tterms = to_dkeys(geo, gamma_interlayer(p), 1)

    for gi in 1:nG
        Gv = bas.G.Qs[gi]
        dk = k .- Gv
        ke = -kin * dot(dk, dk) + p.E0
        for (li, lay) in ((1, :b), (2, :t))
            a = gidx(bas, gi, li)
            H[a, a] += ke + eps / 2 * layer_sign_gamma(lay)
            for (dk_, c) in intraterms
                gj = lookup(bas.G, key_add(geo.Gkey[gi], dk_))
                gj == 0 && continue
                H[a, gidx(bas, gj, li)] += c
            end
        end
        # interlayer: t_Gamma(r) = w_Gamma s_0 is a pure G = 0 transfer
        for (dk_, c) in tterms
            dk_[1] == 0 && dk_[2] == 0 || continue
            H[gidx(bas, gi, 1), gidx(bas, gi, 2)] += c
            H[gidx(bas, gi, 2), gidx(bas, gi, 1)] += conj(c)
        end
    end
    return H
end

"""
    ham_AA_Gamma(geo, bas, p, k; eps = 0.0, convention = :paper)

AA-stacked Gamma-valley Hamiltonian, built from (B32) and (B38):
`V_{Gamma,l}(r)` carries `(-1)^l psi_Gamma` (layer dependent), unlike AB.
Every term is proportional to `s_0`, so each band is two-fold spin degenerate.
"""
function ham_AA_Gamma(geo::MoireGeometry, bas::GammaBasis, p::GammaParams, k;
                      eps::Real = 0.0, convention::Symbol = :paper)
    # layer-dependent phase: the amplitude list depends on (-1)^l
    nG = length(bas.G.Qs)
    H = zeros(ComplexF64, 2nG, 2nG)
    kin = kinetic_prefactor(p.mstar)
    tterms = to_dkeys(geo, gamma_interlayer(p), 1)
    Vby = Dict(l => to_dkeys(geo, aa_gamma_intralayer(geo, p, layer_sign_gamma(l);
                                                       convention = convention), 1)
               for l in (:b, :t))

    for gi in 1:nG
        dk = k .- bas.G.Qs[gi]
        ke = -kin * dot(dk, dk) + p.E0
        for (li, lay) in ((1, :b), (2, :t))
            a = gidx(bas, gi, li)
            H[a, a] += ke + eps / 2 * layer_sign_gamma(lay)
            for (dk_, c) in Vby[lay]
                gj = lookup(bas.G, key_add(geo.Gkey[gi], dk_))
                gj == 0 && continue
                H[a, gidx(bas, gj, li)] += c
            end
        end
        for (dk_, c) in tterms
            dk_[1] == 0 && dk_[2] == 0 || continue
            H[gidx(bas, gi, 1), gidx(bas, gi, 2)] += c
            H[gidx(bas, gi, 2), gidx(bas, gi, 1)] += conj(c)
        end
    end
    return H
end

"""
    ham_AB_Gamma(geo, bas, p, k; eps = 0.0, convention = :paper)

AB-stacked Gamma-valley Hamiltonian of (C11)+(C18) (Eq. (7)).  Here
`V_{Gamma,l}(r)` is layer *independent* (C13, no `(-1)^l`) -- the key difference
from AA -- and `t_Gamma = w_Gamma s_0`.  With the Table II value
`w_Gamma = 300 meV ≫` the kinetic scale, the lowest manifold is a nearly flat
bonding band (see the (C19)-(C25) discussion).
"""
function ham_AB_Gamma(geo::MoireGeometry, bas::GammaBasis, p::GammaParams, k;
                      eps::Real = 0.0, convention::Symbol = :paper)
    intraterms = to_dkeys(geo, ab_gamma_intralayer(geo, p; convention = convention), 1)
    return ham_gamma_spinless(geo, bas, p, intraterms, k; eps = eps)
end

"""promote a spinless Gamma Hamiltonian to the full (layer, spin) space"""
gamma_with_spin(H::AbstractMatrix) = kron(H, Matrix{ComplexF64}(I, 2, 2))

# -----------------------------------------------------------------------------
# unified interface
# -----------------------------------------------------------------------------
"""
    build_basis(geo, model) -> basis object
    hamiltonian(geo, bas, params, model, k; eta, eps, convention) -> Matrix

`model` is one of `:AA_K`, `:AB_K`, `:AA_Gamma`, `:AB_Gamma`.
The returned matrix is Hermitian; `eta` is ignored for the Gamma models (the
Gamma valley is its own manifold).
"""
build_basis(geo::MoireGeometry, model::Symbol) =
    model in (:AA_K, :AB_K) ? kvalley_basis(geo) : gammabasis(geo)

function hamiltonian(geo::MoireGeometry, bas, p, model::Symbol, k;
                     eta::Int = 1, eps::Real = 0.0, convention::Symbol = :paper)
    if model === :AA_K
        return ham_AA_K(geo, bas, p, k; eta = eta, eps = eps, convention = convention)
    elseif model === :AB_K
        return ham_AB_K(geo, bas, p, k; eta = eta, eps = eps, convention = convention)
    elseif model === :AA_Gamma
        return gamma_with_spin(ham_AA_Gamma(geo, bas, p, k; eps = eps, convention = convention))
    elseif model === :AB_Gamma
        return gamma_with_spin(ham_AB_Gamma(geo, bas, p, k; eps = eps, convention = convention))
    else
        error("unknown model $(model)")
    end
end

"""is `H` Hermitian to machine precision?"""
ishermitian_exact(H::AbstractMatrix; tol::Float64 = 1e-10) = norm(H - H') < tol
