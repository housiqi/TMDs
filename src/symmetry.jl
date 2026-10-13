# =============================================================================
# symmetry.jl -- point-group and antiunitary operations as explicit matrices on
#                the plane-wave basis, plus the identity suite that validates
#                the whole calculation.
#
# -----------------------------------------------------------------------------
# How each operator acts on the plane-wave basis
# -----------------------------------------------------------------------------
# Write  c†_{eta,l,r} = sum_p e^{-i p.r} c†_{eta,l,p}.  The paper gives the
# real-space actions (B10), (B18) and (C5); conjugating the expansion is what
# decides whether the momentum label flips:
#
#   C3   c†_{eta,l,p} C3^{-1}   = e^{i eta pi/3} c†_{eta,l, C3 p}       (B10)
#   I    c†_{eta,l,p} I^{-1}    = c†_{eta,lbar, -p}                     (B18)
#   T    c†_{eta,l,p} T^{-1}    = (-eta) c†_{-eta,l,-p}                 (B10)
#   C2yT c†_{eta,l,p} (..)^{-1} = c†_{eta,lbar, C2x p}                  (B10)
#   K    c†_{eta,l,p} K^{-1}    = c†_{eta,lbar, -p}                     (C5)
#
# (I, C2yT and K are antiunitary / unitary exactly as printed; the minus signs
# above are fixed by the conjugation, not chosen.)
#
# The C2yT line deserves a word, because (B10) prints `C2y r` in real space
# while the momentum label above carries `C2x`.  The operation is a *mirror
# composed with time reversal*: on a real-space function f(r),
#
#     (C2y T f)(r) = conj( f(C2y r) ),        C2y : (x,y) -> (-x,y) .
#
# Conjugating the plane-wave expansion e^{-i p.r} turns e^{+i p.r}, so the
# momentum label is reflected twice -- once by C2y and once by the conjugation
# -- i.e. by C2y^{-1} = C2y about the *other* axis.  Explicitly, taking
# C2y r = (-x,y) gives  conj(e^{-i p.(C2y r)}) = conj(e^{-i(p_x(-x)+p_y y)})
# = e^{+i(p_x x - ... )} = e^{-i(C2x p).r}  with  C2x : (x,y) -> (x,-y).
# This is not a convention choice: it is the only reflection that takes
# q1 -> -q1 (q1 is along y), and -q1 lies in the *other* family of (B23), which
# is exactly the layer flip l -> lbar that (B10)/(B14) require.  The reflection
# about y would instead leave q1 (and hence the layer) fixed and is *not* a
# symmetry.  (B4) says the same thing in the monolayer language: in momentum
# space the composite reflects the offset from K, -C2y K_l = K_{-l}, so the two
# K points -- and their layers -- are exchanged.
#
# The basis states are |k,Q> = c†_{eta, l_Q^eta, k-Q} with the plane-wave label
# Q of (B23).  Using l(-Q) = lbar(Q) and l(C2x Q) = lbar(Q) -- both follow from
# (B23), because -q1 and C2x q1 = -q1 both live in the other family -- these
# become
#
#   C3   : |k,Q>_eta -> e^{i eta pi/3} |C3 k,     C3 Q>_eta
#   I    : |k,Q>_eta ->                 | -k,      -Q>_eta
#   T    : |k,Q>_eta -> (-eta)          | -k,      -Q>_{-eta}
#   C2yT : |k,Q>_eta ->                 | C2x k,   C2x Q>_eta
#   K    : |k,Q>_eta ->                 | -k,      -Q>_eta
#
# so that on the matrices h_QQ'(k) of (B27)/(C8)/(C18):
#
#   C3   :  (P3) H(C3 k) (P3)†  = H(k)          P3 : Q -> C3 Q
#   I    :  (Pi) H(-k)  (Pi)†  = H(k)          Pi : Q -> -Q
#   T    :  (Pi) H_{eta}(k) (Pi)† = conj( H_{-eta}(-k) )        (between valleys)
#   C2yT :  (Rx) H(C2x k) (Rx)† = conj( H(k) ) Rx : Q -> C2x Q
#   K    :  (Pi) conj(H(k)) (Pi)† = H(-k)       (AB stacking)
#
# -----------------------------------------------------------------------------
# The Gamma valley: the layer flip has to be put in by hand
# -----------------------------------------------------------------------------
# In the +-K valley the plane-wave label Q = G_M + eta(-1)^l q_1 *encodes* the
# layer, so the label map Q -> -Q carries the layer flip along with it and (B18)
# needs no extra ingredient.  The Gamma-valley label set is just {G_M} -- it
# says nothing about the layer -- so wherever the paper's operation exchanges
# the layers, the permutation must exchange them explicitly:
#
#   (B39)  I  psi†_{r,l} I^{-1} = psi†_{-r,lbar}     (AA Gamma, effective inversion)
#   (C14)  K  psi†_{r,l} K^{-1} = psi†_{r,lbar}      (AB Gamma, effective TR)
#   (B14)  C2yT flips the layer (same as in the +-K valley)
#
# In each case the G_M label goes to -G_M (from r -> -r, or from the
# conjugation in an antiunitary operator), and the *k* argument follows the
# same rule as in the +-K valley: the kinetic term -hbar^2(k-G)^2/2m* is
# invariant under G -> -G only together with k -> -k, so the effective
# inversion / effective TR relate H(k) to H(-k), while the mirror-T composite
# reflects k by C2x (C2y), exactly as above.
#
# Skipping the layer flip is not a harmless omission: at eps = 0 the Gamma
# models are *layer-exchange symmetric on their own* (V is layer independent
# for AB (C13), and for AA at psi_Gamma = 0 the two layers are identical), so
# a permutation without the flip passes trivially -- and then it also passes at
# eps != 0, where the displacement field (-1)^l eps/2 must break both (B40) and
# (C15).  The self-test therefore drives these checks at eps != 0 as well.
#
# Note C3 is *not* a symmetry of H at fixed k: it relates H(C3 k) to H(k).  The
# phase e^{i eta pi/3} is common to all basis states, so it cancels here; it is
# observable only as the eigenvalue of U3 on a C3-invariant eigenstate, which is
# what `c3_rep` checks at Gamma_M.
#
# AB breaking of the effective inversion, and AA of the effective TR, are both
# *reported* rather than asserted -- they are predictions of the paper, and the
# self-test confirms each one comes out the way (B18)/(C5) say it should.
# =============================================================================

# -----------------------------------------------------------------------------
# integer key maps induced by the point-group operations
#
# A key is (n1,n2) with 3Q = n1 b1 + n2 b2; the family (-)^Q = f is read off
# the residue of n1 mod 3 (f = +1  <=>  n1 = 2 mod 3).  Working on keys keeps
# every relabelling exact rather than floating point.
# -----------------------------------------------------------------------------

"""the family (-)^Q = f, read exactly off the key residue"""
@inline key_family(n1::Int) = mod(n1, 3) == 2 ? 1 : -1

"""(m1,m2) of a K-valley key: 3Q = n1 b1 + n2 b2 with Q = m1 b1 + m2 b2 + f q1"""
@inline function key_m(n1::Int, n2::Int)
    f = key_family(n1)
    return (div(n1 + f, 3), div(n2 - 2f, 3), f)
end

"""key of Q = m1 b1 + m2 b2 + f q1"""
@inline key_of(m1::Int, m2::Int, f::Int) = (3m1 - f, 3m2 + 2f)

"""
C3 on a K-valley key.

    C3(m1 b1 + m2 b2 + f q1) = -(m1+m2) b1 + (m1-f) b2 + f q1

because C3 b1 = b2 - b1, C3 b2 = -b1 and C3 q1 = q2 = q1 - b2.  The family, and
hence the layer, is preserved -- as (B14) requires (`h_{eta,l}(C3 r) = h_{eta,l}(r)`).
"""
@inline function c3_key(n1::Int, n2::Int)
    m1, m2, f = key_m(n1, n2)
    return key_of(-(m1 + m2), m1 - f, f)
end

"""inverse of `c3_key`"""
@inline function c3_key_inv(n1::Int, n2::Int)
    m1, m2, f = key_m(n1, n2)
    # invert  (m1,m2) -> (-(m1+m2), m1-f)
    return key_of(m2, -m1 - m2 - f, f)
end

"""inversion: Q -> -Q (the family, i.e. the layer, flips: `l(-Q) = lbar(Q)`)"""
@inline inv_key(n1::Int, n2::Int) = (-n1, -n2)

"""
C2x (the momentum-space form of the mirror in C2yT):  (x,y) -> (x,-y), i.e.

    C2x b1 = b1,  C2x b2 = b1 - b2,  C2x q1 = -q1 ,  so

    C2x(m1 b1 + m2 b2 + f q1) = (m1 + m2) b1 - m2 b2 - f q1
                              = (m1+m2,m2,-f) in the (m1,m2,f) labels
                              = (n1+n2, -n2)   in the integer key.

The family flips, f -> -f, i.e. the layer is exchanged -- as (B14) requires of
C2yT.  (Note C2x b1 = b1 because b1 is along x; C2x b2 = (1/2,-sqrt3/2) = b1-b2.)
"""
@inline function c2x_key(n1::Int, n2::Int)
    m1, m2, f = key_m(n1, n2)
    return key_of(m1 + m2, -m2, -f)
end

# ---- Gamma valley (keys are the (m1,m2) of G_M directly, no family) ---------

"""C3 on a Gamma-valley key: C3 b1 = b2 - b1, C3 b2 = -b1"""
@inline c3_key_g(n1::Int, n2::Int) = (-(n1 + n2), n1)
@inline c3_key_g_inv(n1::Int, n2::Int) = (n2, -n1 - n2)
@inline inv_key_g(n1::Int, n2::Int) = (-n1, -n2)
"""C2x on a Gamma-valley key (no family to flip): b1 -> b1, b2 -> b1 - b2"""
@inline c2x_key_g(n1::Int, n2::Int) = (n1 + n2, -n2)

"""
C2y on a Gamma-valley key: C2y b1 = -b1, C2y b2 = b2 - b1 (the mirror about y
in the paper's axis convention, as used in (B14)/(B40)).

C2yT and C2xT are two *different* composites, and the G_M label map is not the
same for the two: taking (T M f)(r) = conj(f(M r)) with f = exp(+i(k-G).r) gives
exp(-i M(k-G).r) = exp(+i(M G - M k).r), so

    C2yT :  G -> -C2y G = C2x G   and   k -> -C2y k = C2x k     (c2x_key_g)
    C2xT :  G -> -C2x G = C2y G   and   k -> -C2x k = C2y k     (c2y_key_g)

(using -C2y = C2x as 2x2 matrices).  AA uses C2yT (B14) and AB uses C2xT; the
two maps are each other's inverse, so mixing them up would flip the mirror
axis.  They are distinct precisely when the model violates C2x (e.g. AA at
psi != 0), which is why the self-test checks each model against its own.
"""
@inline c2y_key_g(n1::Int, n2::Int) = (-n1 - n2, n2)

# -----------------------------------------------------------------------------
# permutation matrices
# -----------------------------------------------------------------------------

"""the key list of a Gamma basis, indexed as the basis is"""
function keylist(bas::GammaValleyBasis)
    ks = Vector{NTuple{2,Int}}(undef, length(bas.Qs))
    for (k, i) in bas.index
        ks[i] = k
    end
    return ks
end

"""permutation array `p` with `p[a]` = index of `keymap(Q_a)`; errors if it leaves the set"""
function basis_perm(bas::KValleyBasis, keymap::Function)
    p = Vector{Int}(undef, length(bas.Qs))
    for a in eachindex(bas.Qs)
        b = lookup(bas, keymap(bas.key[a]...))
        b == 0 && error("symmetry leaves the K-valley basis at key $(bas.key[a])")
        p[a] = b
    end
    return p
end

function basis_perm(bas::GammaValleyBasis, keymap::Function)
    ks = keylist(bas)
    p = Vector{Int}(undef, length(bas.Qs))
    for a in eachindex(bas.Qs)
        b = get(bas.index, keymap(ks[a]...), 0)
        b == 0 && error("symmetry leaves the Gamma basis at key $(ks[a])")
        p[a] = b
    end
    return p
end

"""same, promoted to the full (plane wave, layer, spin) space of the Gamma models

`layerswap = true` exchanges b <-> t on the way.  It is *required* whenever the
paper's operation flips the layer -- (B39) for the AA effective inversion, (C14)
for the AB effective TR, (B14)/(C15) for the mirror-T composites -- because the
Gamma-valley plane-wave label carries no layer information (see the header).
The spin index is left alone: (B39)/(C14)/C3 all leave s alone (C9).
"""
function basis_perm(bas::GammaBasis, keymap::Function; layerswap::Bool = false)
    pG = basis_perm(bas.G, keymap)
    nG = length(bas.G.Qs)
    p = Vector{Int}(undef, 4nG)
    for gi in 1:nG, li in 1:2, si in 1:2
        lj = layerswap ? 3 - li : li          # li = 1 <-> b, 2 <-> t
        p[((gi - 1) * 2 + (li - 1)) * 2 + si] = ((pG[gi] - 1) * 2 + (lj - 1)) * 2 + si
    end
    return p
end

"""permute a matrix: `H[p,p]`"""
permute(H::AbstractMatrix, p::Vector{Int}) = H[p, p]

"""C3 permutation on the *spinless* (G_M, layer) basis of the Gamma models

C3 rotates r, so it does not exchange the layers -- no `layerswap` here, unlike
(B39)/(C14)/(B14).  This is the basis in which the paper quotes the spinless C3
eigenvalues of (C24).
"""
function c3_perm_spinless(bas::GammaBasis)
    pG = basis_perm(bas.G, c3_key_g)
    nG = length(bas.G.Qs)
    p = Vector{Int}(undef, 2nG)
    for gi in 1:nG, li in 1:2
        p[gidx(bas, gi, li)] = gidx(bas, pG[gi], li)
    end
    return p
end

"""
    top_multiplets(H, P; nlev = 2, etol = 1e-6) -> Vector{@NamedTuple{E,mult,c3}}

The `nlev` highest distinct levels of `H`, each as its energy, its multiplicity
(number of states within `etol`) and the eigenvalues of the operator `P`
restricted to that degenerate subspace.  This is how the paper's (C24)
predictions -- multiplicities and spinless C3 eigenvalues of the Gamma-valley
flat bands -- are turned into numbers that can be checked.
"""
function top_multiplets(H::AbstractMatrix, P::AbstractMatrix; nlev::Int = 2,
                        etol::Float64 = 1e-6)
    E, V = eigen(Hermitian(Matrix(H)))
    out = NamedTuple{(:E, :mult, :c3),Tuple{Float64,Int,Vector{ComplexF64}}}[]
    i = length(E)
    while i >= 1 && length(out) < nlev
        j = i
        while j > 1 && abs(E[j] - E[j-1]) < etol
            j -= 1
        end
        push!(out, (E = E[i], mult = i - j + 1,
                    c3 = eigvals(V[:, j:i]' * P * V[:, j:i])))
        i = j - 1
    end
    return out
end

# -----------------------------------------------------------------------------
# the checks
# -----------------------------------------------------------------------------
"""one line of the self-test table"""
struct CheckResult
    name::String
    value::Float64
    passed::Bool
    expected::Bool          # whether this identity is supposed to hold
    note::String
end

"""
    c3_rep(H, P3) -> Vector{ComplexF64}

Eigenvalues of the C3 permutation restricted to each degenerate eigenspace of
`H`.  For a C3-invariant k-point (Gamma_M) every one of them must be a cube
root of unity -- that is the C3 representation content of the bands.
"""
function c3_rep(H::AbstractMatrix, P3::AbstractMatrix; etol::Float64 = 1e-7)
    E, V = eigen(Hermitian(Matrix(H)))
    n = size(H, 1)
    reps = ComplexF64[]
    i = 1
    while i <= n
        j = i
        while j < n && abs(E[j + 1] - E[i]) < etol
            j += 1
        end
        B = V[:, i:j]' * P3 * V[:, i:j]
        append!(reps, eigvals(B))
        i = j + 1
    end
    return reps
end

"""
    pair_gap(E) -> Float64

Largest `E[2i] - E[2i-1]` for an ascending spectrum `E`.  Zero iff the spectrum
is exactly two-fold degenerate throughout -- the natural way to *measure* a
degeneracy the paper predicts, instead of counting tolerance windows.
"""
function pair_gap(E::Vector{Float64})
    m = 0.0
    for i in 1:2:(length(E) - 1)
        m = max(m, E[i+1] - E[i])
    end
    return m
end

"""how many eigenvalues of `E` (strictly ascending) are two-fold degenerate"""
function pair_degeneracy(E::Vector{Float64}; tol::Float64 = 1e-6)
    n = length(E)
    npair = 0
    i = 1
    while i < n
        if abs(E[i + 1] - E[i]) < tol
            npair += 1
            i += 2
        else
            i += 1
        end
    end
    return npair
end

"""
    symmetry_report(geo; eps = 0.0, Nshell = 5) -> Vector{CheckResult}

Run the whole identity suite.  `expected = false` entries are identities the
paper says should *fail* (e.g. effective inversion for AB stacking); they are
reported so the failure is visible evidence rather than silence.
"""
function symmetry_report(geo::MoireGeometry; eps::Real = 0.0, kgen = [0.0173, -0.0091])
    out = CheckResult[]
    tol = 1e-9
    k = collect(float.(kgen))
    c3k = C3M * k
    C2xM = [1.0 0.0; 0.0 -1.0]
    c2xk = C2xM * k

    specs = [(m, p) for (m, p) in
             ((:AA_K, AA_K_FH), (:AA_K, AA_K_FHSH), (:AB_K, AB_K),
              (:AA_Gamma, AA_GAMMA), (:AB_Gamma, AB_GAMMA))]
    # AA_K appears twice; label the two parameter sets apart
    labels = Dict(:AA_K => "AA±K", :AB_K => "AB±K",
                  :AA_Gamma => "AA-Γ", :AB_Gamma => "AB-Γ")

    for (model, prm) in specs
        ptag = model === :AA_K ? " ($(label(prm)))" : ""
        bas = build_basis(geo, model)
        tag = labels[model] * ptag
        ketas = model in (:AA_K, :AB_K) ? (1, -1) : (1,)

        # --- 1. Hermiticity -------------------------------------------------
        let worst = 0.0
            for kk in (k, c3k), eta in ketas
                H = hamiltonian(geo, bas, prm, model, kk; eta = eta, eps = eps)
                worst = max(worst, norm(H - H'))
            end
            push!(out, CheckResult("$tag: Hermiticity", worst, worst < tol, true,
                                   "‖H - H'‖"))
        end

        # --- 2. :paper vs :transposed (plain transpose => identical spectra) --
        # (B22) is a field-operator expansion, so the state carries the
        # conjugate kernel and the rule is h_{Q,Q'} = O_{eta, eta(Q-Q')} -- the
        # printed (B27)/(C8)/(C18).  Reading the kernel as the wavefunction
        # instead gives h_{Q,Q'} = O_{eta, eta(Q'-Q)}, the plain transpose.  For
        # Hermitian H the two differ by psi -> -psi and have identical spectra,
        # so this is an internal consistency check of the amplitude assignment.
        let worst = 0.0, same = true
            for eta in ketas
                Hp = hamiltonian(geo, bas, prm, model, k; eta = eta, eps = eps,
                                 convention = :paper)
                Hr = hamiltonian(geo, bas, prm, model, k; eta = eta, eps = eps,
                                 convention = :transposed)
                worst = max(worst, norm(Hp - transpose(Hr)))
                same &= norm(Hp - Hr) < tol
            end
            push!(out, CheckResult("$tag: H_paper == transpose(H_transposed)", worst,
                                   worst < tol, true,
                                   same ? "psi = 0: conventions coincide, check degenerate" :
                                          "psi != 0: conventions genuinely differ"))
        end

        # --- 3. C3 covariance ----------------------------------------------
        P3 = basis_perm(bas, model === :AA_K || model === :AB_K ? c3_key : c3_key_g)
        P3m = Matrix{Float64}(I, length(P3), length(P3))[P3, :]   # P3 H(C3k) P3' = H(k)
        let worst = 0.0
            for eta in ketas
                Hk = hamiltonian(geo, bas, prm, model, k; eta = eta, eps = eps)
                Hc = hamiltonian(geo, bas, prm, model, c3k; eta = eta, eps = eps)
                worst = max(worst, norm(permute(Hc, P3) - Hk))
            end
            push!(out, CheckResult("$tag: H(C3k) = P3 H(k) P3'", worst, worst < tol,
                                   true, "C3 covariance at generic k"))
        end

        # --- 4. C3 representation at Gamma_M --------------------------------
        let worst = 0.0
            for eta in ketas
                H0 = hamiltonian(geo, bas, prm, model, [0.0, 0.0]; eta = eta, eps = eps)
                reps = c3_rep(H0, P3m)
                worst = max(worst, maximum(abs.(reps .^ 3 .- 1)))
                worst = max(worst, maximum(abs.(abs.(reps) .- 1)))
            end
            push!(out, CheckResult("$tag: C3 reps at Γ_M ∈ {1,ω,ω²}", worst,
                                   worst < 1e-6, true, "e^{iηπ/3} factors out at fixed k"))
        end

        # --- 5/6. effective inversion ---------------------------------------
        # (B18) in the +-K valley, (B39) in the AA Gamma valley:
        # I flips the layer AND the sign of r.  The +-K-valley label
        # Q -> -Q flips the family, hence the layer, by itself; the Gamma-valley
        # label set carries no layer at all, so the swap is explicit there.
        Pi = model in (:AA_K, :AB_K) ? basis_perm(bas, inv_key) :
                                       basis_perm(bas, inv_key_g; layerswap = true)
        let worst = 0.0
            for eta in ketas
                Hk = hamiltonian(geo, bas, prm, model, k; eta = eta, eps = eps)
                Hm = hamiltonian(geo, bas, prm, model, -k; eta = eta, eps = eps)
                worst = max(worst, norm(permute(Hm, Pi) - Hk))
            end
            should = model in (:AA_K, :AA_Gamma, :AB_Gamma) && eps == 0
            push!(out, CheckResult("$tag: effective inversion P H(-k) P' = H(k)",
                                   worst, worst < tol, should,
                                   model === :AB_Gamma ?
                                       "(C13) V is layer independent => the model is " *
                                       "layer-exchange symmetric; with eps = 0 an " *
                                       "effective inversion follows for free" :
                                   (should ? "(B18)/(B39): flips layer and r" :
                                             "(B19)/(C8)/(C15): broken — expected")))
        end

        # --- 7. mirror x T ---------------------------------------------------
        # The real-space generator is the mirror C2y (AA) or C2x (AB), but the
        # composite mirror-T acts on the plane-wave label as the *other*
        # reflection: the conjugation in T reflects the label a second time
        # (see the header and the c2y_key_g docstring).  Both composites flip
        # the layer, which the Gamma-valley permutation has to do by hand.
        if model === :AA_K
            let worst = 0.0, Rx = basis_perm(bas, c2x_key)
                for eta in ketas
                    Hk = hamiltonian(geo, bas, prm, model, k; eta = eta, eps = eps)
                    Hc = hamiltonian(geo, bas, prm, model, c2xk; eta = eta, eps = eps)
                    worst = max(worst, norm(permute(Hc, Rx) - conj.(Hk)))
                end
                push!(out, CheckResult("$tag: C2yT: Rx H(C2x k) Rx' = conj H(k)",
                                       worst, worst < tol, eps == 0,
                                       "(B10)/(B14); must break at eps != 0"))
            end
        elseif model === :AA_Gamma
            let worst = 0.0, Rx = basis_perm(bas, c2x_key_g; layerswap = true)
                Hk = hamiltonian(geo, bas, prm, model, k; eps = eps)
                Hc = hamiltonian(geo, bas, prm, model, c2xk; eps = eps)
                worst = norm(permute(Hc, Rx) - conj.(Hk))
                push!(out, CheckResult("$tag: C2yT: Rx H(C2x k) Rx' = conj H(k)",
                                       worst, worst < tol, eps == 0,
                                       "(B14); must break at eps != 0"))
            end
        elseif model === :AB_Gamma
            # AB uses the other mirror, C2x (text after (C15)): the composite
            # C2xT reflects the label by C2y and k by C2y as well.
            let worst = 0.0, Rx = basis_perm(bas, c2y_key_g; layerswap = true)
                Hk = hamiltonian(geo, bas, prm, model, k; eps = eps)
                Hc = hamiltonian(geo, bas, prm, model, [-k[1], k[2]]; eps = eps)
                worst = norm(permute(Hc, Rx) - conj.(Hk))
                push!(out, CheckResult("$tag: C2xT: Rx H(C2y k) Rx' = conj H(k)",
                                       worst, worst < tol, eps == 0,
                                       "(C15); must break at eps != 0"))
            end
        end

        # --- 8. valley TR ----------------------------------------------------
        if model in (:AA_K, :AB_K)
            let worst = 0.0
                Hp = hamiltonian(geo, bas, prm, model, -k; eta = 1, eps = eps)
                Hm = hamiltonian(geo, bas, prm, model, k; eta = -1, eps = eps)
                worst = norm(permute(Hm, Pi) - conj.(Hp))
                push!(out, CheckResult("$tag: T: Pi H_{-η}(k) Pi' = conj H_η(-k)",
                                       worst, worst < tol, true, "(B10)"))
            end

            # --- 9. AB effective TR K (C5) ----------------------------------
            let worst = 0.0
                Hk = hamiltonian(geo, bas, prm, model, k; eta = 1, eps = eps)
                Hm = hamiltonian(geo, bas, prm, model, -k; eta = 1, eps = eps)
                worst = norm(permute(Hm, Pi) - conj.(Hk))
                should = model === :AB_K && eps == 0
                push!(out, CheckResult("$tag: effective K: Pi conj H(k) Pi' = H(-k)",
                                       worst, worst < tol, should,
                                       should ? "(C5) AB lowest harmonics" :
                                                "(B18)/(C8): not a symmetry of AA"))
            end

            # --- 10. the two valleys -----------------------------------------
            let Hp = hamiltonian(geo, bas, prm, model, k; eta = 1, eps = eps),
                Hm = hamiltonian(geo, bas, prm, model, k; eta = -1, eps = eps)

                if model === :AB_K && eps == 0
                    v = norm(Hp - Hm)
                    push!(out, CheckResult("$tag: H_{+1}(k) == H_{-1}(k) identically",
                                           v, v < tol, true,
                                           "AB ±K has no η-dependence at ε = 0"))
                else
                    v = maximum(abs.(sort(real.(eigvals(Hermitian(Matrix(Hp))))) .-
                                      sort(real.(eigvals(Hermitian(Matrix(Hm)))))))
                    push!(out, CheckResult("$tag: spec H_{+1}(k) == spec H_{-1}(k)",
                                           v, v < 1e-7, true,
                                           "AA: via effective inversion + ψ → -ψ"))
                end
            end
        end

        # --- 10b. AB Gamma effective TR (C14) --------------------------------
        # K psi†_{r,l} K^{-1} = psi†_{r,lbar}: antiunitary, flips the layer and
        # (through the conjugation and the kinetic term) maps k -> -k.  This is
        # the Gamma-valley analogue of the +-K-valley (C5), and it too needs the
        # explicit layer swap.
        if model === :AB_Gamma
            let worst = 0.0, Pk = basis_perm(bas, inv_key_g; layerswap = true)
                Hk = hamiltonian(geo, bas, prm, model, k; eps = eps)
                Hm = hamiltonian(geo, bas, prm, model, -k; eps = eps)
                worst = norm(permute(Hm, Pk) - conj.(Hk))
                push!(out, CheckResult("$tag: effective K (C14): P conj H(k) P' = H(-k)",
                                       worst, worst < tol, eps == 0,
                                       "(C14)/(C15): broken by eps — expected"))
            end
        end

        # --- 11a. tunneling is strictly interlayer (B13) ---------------------
        # (B13) gives h_{eta,l}(r) Hermitian and *layer diagonal*, so every
        # off-diagonal-in-layer entry must come from t_eta(r), whose transfers
        # are +-q_i and +-q_{2i}.  A same-family entry at one of those transfers
        # would be an error.  (The honeycomb structure already forbids it: Q
        # +- q_i either lands on the other sublattice or off the honeycomb
        # entirely.  Checking it also guards the amplitude bookkeeping, which is
        # what the earlier probes of this file were worried about.)
        if model in (:AA_K, :AB_K)
            tdk = to_dkeys(geo, model === :AA_K ? aa_k_interlayer(geo, prm, 1) :
                                             ab_k_interlayer(geo, prm), 3)
            tdkd = [d for (d, c) in tdk]
            let worst = 0.0
                for eta in ketas
                    H = hamiltonian(geo, bas, prm, model, k; eta = eta, eps = eps)
                    for a in eachindex(bas.f), b in eachindex(bas.f)
                        (a == b || bas.f[a] != bas.f[b]) && continue
                        d = (bas.key[b][1] - bas.key[a][1], bas.key[b][2] - bas.key[a][2])
                        any(dd -> dd == d, tdkd) && (worst = max(worst, abs(H[a, b])))
                    end
                end
                push!(out, CheckResult("$tag: tunneling strictly interlayer (B13)",
                                       worst, worst < tol, true,
                                       "no same-layer ±q_i / ±q_2i entries"))
            end
        end

        # --- 11b. AB +-K: the two layer blocks (w = 0) ------------------------
        # Table II sets w = 0, so the layers decouple: H is exactly block
        # diagonal in the family index, and the two blocks are related by the
        # mirror, spec H_t(k) = spec H_b(C2x k).  A layer degeneracy is
        # therefore exact exactly when C2x k lies in the same C3 orbit as k --
        # on the 60-degree ray Gamma_M -> M_M, where C2x = C3^2 -- and absent
        # off it, e.g. along Gamma_M -> K_M.  That is the paper's Fig. 4(c)
        # pattern: fourfold (layer x spin) along Gamma_M-M_M, twofold (spin
        # only) at K_M.  (An earlier version of this check tested a generic k
        # and was ill-posed: the blocks coincide only on the mirror line.)
        if model === :AB_K
            fb = findall(==(1), bas.f)
            ft = findall(==(-1), bas.f)
            spec(M) = sort(real.(eigvals(Hermitian(Matrix(M)))))
            Hb(kk) = hamiltonian(geo, bas, prm, model, kk; eps = 0.0)[fb, fb]
            Ht(kk) = hamiltonian(geo, bas, prm, model, kk; eps = 0.0)[ft, ft]
            H0 = hamiltonian(geo, bas, prm, model, k; eps = 0.0)
            let offblk = maximum([abs(H0[a, b]) for a in fb, b in ft])
                push!(out, CheckResult("$tag: H exactly block diagonal in the layer (w=0)",
                                       offblk, offblk < tol, true,
                                       "the two layers decouple"))
            end
            let v = maximum(abs.(spec(Ht(k)) .- spec(Hb([k[1], -k[2]]))))
                push!(out, CheckResult("$tag: spec H_t(k) == spec H_b(C2x k) (w=0)",
                                       v, v < tol, true, "mirror-related layer blocks"))
            end
            segdeg(a, b, n) = [pair_gap(diagonalize_point(geo, bas, prm, model,
                                                          a .+ (b .- a) .* t))
                               for t in range(0.0, 1.0; length = n)]
            # Gamma_M -> K_M: skip the shared endpoint Gamma_M, where the
            # degeneracy is exact for the trivial reason k = C2x k = 0
            onGM = [pair_gap(diagonalize_point(geo, bas, prm, model, geo.K_M .* t))
                    for t in range(0.2, 1.0; length = 5)]
            onMG = segdeg(geo.Mp_M, [0.0, 0.0], 6)     # M_M -> Gamma_M
            push!(out, CheckResult("$tag: layer degeneracy on Γ_M–M_M (w=0)",
                                   maximum(onMG), maximum(onMG) < 1e-6, true,
                                   "C2x = C3² along the 60° ray"))
            push!(out, CheckResult("$tag: no layer degeneracy on Γ_M–K_M (w=0)",
                                   minimum(onGM), minimum(onGM) > 1e-3, false,
                                   "C2x k falls in the other K class there"))
        end
        if model in (:AA_Gamma, :AB_Gamma)
            H = hamiltonian(geo, bas, prm, model, k; eps = eps)
            E = sort(real.(eigvals(Hermitian(Matrix(H)))))
            np = pair_degeneracy(E)
            ok = np == length(E) ÷ 2
            push!(out, CheckResult("$tag: exact 2-fold spin degeneracy",
                                   Float64(np), ok, true,
                                   "every term ∝ s_0, $np of $(length(E) ÷ 2) pairs"))

            # --- 12. the (C24) harmonic-oscillator multiplet at Gamma_M --------
            # With m_Gamma* = 10 m_e the moire potential -V_Gamma(r) confines the
            # top valence states to |r| ~ 3.7 Ang, far inside the 52 Ang moire
            # cell, so they are the (n_x,n_y) eigenstates of an array of
            # decoupled harmonic oscillators (App. C 3).  (C24) then predicts
            # that the highest valence level is (0,0) -- one state spinless,
            # spinless C3 eigenvalue 1 -- and the second is the (1,0),(0,1)
            # doublet with C3 = e^{+-i2pi/3}; including spin, the paper's 2-fold
            # and 4-fold degeneracies (Fig. 4(c) caption, and (i)/(ii) after
            # (C25)).  The manifolds must also be ultra-flat over the whole MBZ:
            # that localization is exactly what "ultra-flat" means here.
            spinless = model === :AB_Gamma ?
                       ham_AB_Gamma(geo, bas, prm, [0.0, 0.0]; eps = eps) :
                       ham_AA_Gamma(geo, bas, prm, [0.0, 0.0]; eps = eps)
            Pc3 = Matrix{Float64}(I, 2 * length(bas.G.Qs), 2 * length(bas.G.Qs))[
                      c3_perm_spinless(bas), :]
            lev = top_multiplets(spinless, Pc3)
            l0, l1 = lev[1], lev[2]
            push!(out, CheckResult("$tag: (C24) top level: 1 state, spinless C3 = 1",
                                   max(abs(l0.c3[1] - 1), Float64(abs(l0.mult - 1))),
                                   l0.mult == 1 && abs(l0.c3[1] - 1) < 1e-8, true,
                                   "(n_x,n_y) = (0,0); 2-fold once spin is included"))
            push!(out, CheckResult("$tag: (C24) next level: 2 states, C3 = e^{±i2π/3}",
                                   max(maximum(abs.(l1.c3 .^ 3 .- 1)),
                                       Float64(abs(l1.mult - 2))),
                                   l1.mult == 2 &&
                                   maximum(abs.(l1.c3 .^ 3 .- 1)) < 1e-7 &&
                                   minimum(abs.(l1.c3 .- 1)) > 0.5, true,
                                   "(1,0),(0,1); 4-fold once spin is included"))
            let hw = hbar_omega(geo, prm), gap = l0.E - l1.E
                push!(out, CheckResult("$tag: (C24)/(C25) flat-band splitting ≈ ħω",
                                       abs(gap - hw), abs(gap - hw) < 3.0, true,
                                       "ħω = $(round(hw, digits = 2)) meV; the exact " *
                                       "lattice value is $(round(gap, digits = 2)) meV — " *
                                       "smaller by the quartic anharmonicity of " *
                                       "-V_Γ(r), whose first-order perturbation " *
                                       "accounts for it (see README/App. C 3)"))
            end
            if model === :AB_Gamma
                # E_Gamma is not tabulated (App. B 3); the default in
                # parameters.jl is calibrated so that the flat band sits 30 meV
                # under the +-K valence band maximum (Sec. III.C, Fig. 4(c)).
                let path = bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]];
                                      nseg = [12, 6, 12]),
                    pos = vbm(geo, bas, prm, :AB_Gamma, path) -
                          vbm(geo, build_basis(geo, :AB_K), AB_K, :AB_K, path)
                    push!(out, CheckResult("$tag: flat band sits 30 meV below the ±K VBM",
                                           abs(pos + 30), abs(pos + 30) < 8, true,
                                           "measured $(round(pos, digits = 2)) meV " *
                                           "(sum along the path; E_Γ is a free parameter)"))
                end
            end
        end
    end

    # --- 13. the harmonic scale of the AB Gamma valley -----------------------
    let w = hbar_omega(geo, AB_GAMMA)
        push!(out, CheckResult("AB-Γ: ħω from (C25) [meV]", w, abs(w - 56.76) < 0.05,
                               true, "paper quotes 56.76 meV"))
    end

    return out
end

"""
    print_report(results; io = stdout)

Pretty-print the self-test table; returns `true` if every `expected` identity
holds.
"""
function print_report(results::Vector{CheckResult}; io::IO = stdout)
    println(io, "="^78)
    println(io, "  symmetry / consistency self-test    (θ = $(THETA_DEFAULT_DEG)°, ε = 0)")
    println(io, "="^78)
    @printf(io, "  %-46s %10s  %s\n", "check", "residual", "status")
    println(io, "-"^78)
    ok = true
    for c in results
        if c.expected
            @printf(io, "  %-46s %10.2e  %s\n", c.name, c.value, c.passed ? "ok" : "FAIL")
            ok &= c.passed
        else
            @printf(io, "  %-46s %10.2e  %s\n", c.name, c.value,
                    c.passed ? "flag" : "n/a")
            println(io, "      ↳ not a symmetry of this model — " * c.note)
        end
    end
    println(io, "="^78)
    println(io, ok ? "  ALL EXPECTED IDENTITIES HOLD" : "  SOME EXPECTED IDENTITIES FAILED")
    println(io, "="^78)
    return ok
end
