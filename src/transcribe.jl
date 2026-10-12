# =============================================================================
# transcribe.jl -- verbatim transcription of the paper's *momentum-space*
#                  equations, as an independent cross-check of kmodels.jl
#
# This file deliberately does NOT reuse the integer-key machinery of
# geometry.jl or the Fourier-component tables of potentials.jl.  It scans every
# pair of basis states and tests the printed delta functions directly on the
# float Q vectors.  Agreement with kmodels.jl therefore validates both the key
# bookkeeping and the Fourier-amplitude bookkeeping.
#
# Equations transcribed:
#   (B27)  AA stacking, +-K valley
#   (C8)   AB stacking, +-K valley
#   (C18)  AB stacking, Gamma valley
#
# There is no printed *momentum-space* form for the AA Gamma valley (App. B 3
# only gives the real-space (B32)/(B38)); `transcribe_hamiltonian` refuses that
# model rather than inventing one.
# =============================================================================

const DTOL = 1e-8

"""float delta: are the two 2-vectors equal within DTOL?"""
same_vec(u, v; tol::Float64 = DTOL) = norm(u .- v) < tol

"""
    transcribe_hamiltonian(geo, bas, p, model, k; eta = 1, eps = 0.0)

Direct transcription of the printed momentum-space Hamiltonian.  Returns a
matrix with the *same* shape as `hamiltonian(...)` for the same arguments, so
that the two can be compared element by element: the ``+/-K`` models have no
spin index, while the Gamma models are promoted with `gamma_with_spin`.
"""
function transcribe_hamiltonian(geo::MoireGeometry, bas, p, model::Symbol, k;
                                eta::Int = 1, eps::Real = 0.0)
    if model === :AA_K
        return transcribe_AA_K(geo, bas, p, k; eta = eta, eps = eps)
    elseif model === :AB_K
        return transcribe_AB_K(geo, bas, p, k; eta = eta, eps = eps)
    elseif model === :AB_Gamma
        return gamma_with_spin(transcribe_AB_Gamma(geo, bas, p, k; eps = eps))
    elseif model === :AA_Gamma
        error("App. B 3 gives no momentum-space Gamma form for AA stacking; " *
              "use ham_AA_Gamma (built from (B32)+(B38)) instead.")
    else
        error("unknown model $(model)")
    end
end

# --- (B27) ------------------------------------------------------------------
"""
    transcribe_AA_K(geo, bas, p, k; eta, eps)

Eq. (B27), term by term:

    h_{Q Q'} = delta_QQ' [ -hbar^2 (k-Q)^2/(2m*) + (eps/2) eta (-)^Q ]
             + V  sum_i [ exp(-eta (-)^Q i psi) delta_{Q+g_i,   Q'}
                        + exp(+eta (-)^Q i psi) delta_{Q-g_i,   Q'} ]
             + V2 sum_i [ delta_{Q+g_{2i},Q'} + delta_{Q,Q'+g_{2i}} ]
             + w  sum_i [ delta_{Q+q_i,  Q'} + delta_{Q,Q'+q_i}  ]
             + w2 sum_i [ delta_{Q+q_{2i},Q'} + delta_{Q,Q'+q_{2i}} ]

with `(-)^Q = f`, the Q family (the paper's `(-)^Q`, fixed by `eta (-)^Q =
(-1)^l` in (B26)).
"""
function transcribe_AA_K(geo::MoireGeometry, bas::KValleyBasis, p::KParams, k;
                         eta::Int = 1, eps::Real = 0.0)
    n = length(bas.Qs)
    H = zeros(ComplexF64, n, n)
    kin = kinetic_prefactor(p.mstar)
    for a in 1:n, b in 1:n
        Q, Qp = bas.Qs[a], bas.Qs[b]
        f = bas.f[a]                       # the paper's (-)^Q
        if a == b
            dk = k .- Q
            H[a, b] += -kin * dot(dk, dk) + eps / 2 * eta * f
        end
        for i in 1:3
            same_vec(Q .+ geo.g[i], Qp) && (H[a, b] += p.V * exp(-im * eta * f * p.psi))
            same_vec(Q .- geo.g[i], Qp) && (H[a, b] += p.V * exp(+im * eta * f * p.psi))
            (same_vec(Q .+ geo.g2[i], Qp) || same_vec(Q, Qp .+ geo.g2[i])) && (H[a, b] += p.V2)
            (same_vec(Q .+ geo.q[i], Qp) || same_vec(Q, Qp .+ geo.q[i])) && (H[a, b] += p.w)
            (same_vec(Q .+ geo.q2[i], Qp) || same_vec(Q, Qp .+ geo.q2[i])) && (H[a, b] += p.w2)
        end
    end
    return H
end

# --- (C8) -------------------------------------------------------------------
"""
    transcribe_AB_K(geo, bas, p, k; eta, eps)

Eq. (C8):

    h_{Q Q'} = delta_QQ' [ -hbar^2 (k-Q)^2/(2m*) + (eps/2) eta (-)^Q ]
             + V  sum_i sum_{alpha=+-1} exp(-alpha i psi)      delta_{Q+alpha g_i, Q'}
             + w  sum_i sum_{alpha=+-1} exp(i alpha (i-1) 2pi/3) delta_{Q, Q'+alpha q_i}
"""
function transcribe_AB_K(geo::MoireGeometry, bas::KValleyBasis, p::KParams, k;
                         eta::Int = 1, eps::Real = 0.0)
    n = length(bas.Qs)
    H = zeros(ComplexF64, n, n)
    kin = kinetic_prefactor(p.mstar)
    for a in 1:n, b in 1:n
        Q, Qp = bas.Qs[a], bas.Qs[b]
        f = bas.f[a]
        if a == b
            dk = k .- Q
            H[a, b] += -kin * dot(dk, dk) + eps / 2 * eta * f
        end
        for i in 1:3, alpha in (1, -1)
            same_vec(Q .+ alpha .* geo.g[i], Qp) && (H[a, b] += p.V * exp(-im * alpha * p.psi))
            same_vec(Q, Qp .+ alpha .* geo.q[i]) && (H[a, b] += p.w * exp(im * alpha * (i - 1) * 2pi / 3))
        end
    end
    return H
end

# --- (C18) ------------------------------------------------------------------
"""
    transcribe_AB_Gamma(geo, bas, p, k; eps)

Eq. (C18), per layer:

    h_{ll, G G'} = [ -hbar^2 (k-G)^2/(2m_Gamma*) + (-)^l eps/2 + E_Gamma ] delta_{G,G'}
                 + V_Gamma sum_i [ exp(+i psi_Gamma) delta_{G, G'-g_i}
                                 + exp(-i psi_Gamma) delta_{G, G'+g_i} ]
    h_{l lbar, G G'} = w_Gamma delta_{G,G'}

(The last line is printed in (C18) with the same `h_{ll,...}` subscript as the
previous one -- an obvious typesetting slip; `l lbar` is what the tunnel term
must be.)  Returns the spinless 2-layer matrix; promote with `gamma_with_spin`.
"""
function transcribe_AB_Gamma(geo::MoireGeometry, bas::GammaBasis, p::GammaParams, k;
                             eps::Real = 0.0)
    nG = length(bas.G.Qs)
    H = zeros(ComplexF64, 2nG, 2nG)
    kin = kinetic_prefactor(p.mstar)
    for gi in 1:nG, gj in 1:nG
        G, Gp = bas.G.Qs[gi], bas.G.Qs[gj]
        for (li, lay) in ((1, :b), (2, :t))
            a, b = gidx(bas, gi, li), gidx(bas, gj, li)
            if gi == gj
                dk = k .- G
                H[a, b] += -kin * dot(dk, dk) + layer_sign_gamma(lay) * eps / 2 + p.E0
            end
            for i in 1:3
                same_vec(G, Gp .- geo.g[i]) && (H[a, b] += p.V * exp(+im * p.psi))
                same_vec(G, Gp .+ geo.g[i]) && (H[a, b] += p.V * exp(-im * p.psi))
            end
        end
        # interlayer
        if gi == gj
            H[gidx(bas, gi, 1), gidx(bas, gj, 2)] += p.w
            H[gidx(bas, gi, 2), gidx(bas, gj, 1)] += p.w
        end
    end
    return H
end
