# =============================================================================
# geometry.jl -- moire geometry of twisted bilayer TMDs
# Every formula below cites PRB 109, 205121 (2024) by equation number.
# Lengths in Angstrom, momenta in Angstrom^-1, energies in meV.
#   (B1)   a_M   = a_0 / (2 sin(theta/2))
#   (2)*   b_{M,1} = 4pi/(sqrt(3) a_M) (1,0)^T
#   (2)*   b_{M,2} = 4pi/(sqrt(3) a_M) (1/2, sqrt(3)/2)^T
#   (2)*   a_{M,1} = a_M (sqrt(3)/2, -1/2)^T   ,  a_{M,2} = C3 a_{M,1}
#          g_i   = C3^{i-1} b_{M,1}                        (i = 1,2,3)
#          g_{2,1} = b_{M,1} + b_{M,2} ,  g_{2,i} = C3^{i-1} g_{2,1}
#   (B12)  q_1 = K_b - K_t = 4pi/(3a_0) 2 sin(theta/2) (0,1)^T
#          q_i   = C3^{i-1} q_1
#          q_{2,1} = b_{M,1} + q_1 ,  q_{2,i} = C3^{i-1} q_{2,1}
#   (B23)  Q_l^eta = { G_M + eta (-1)^l q_1 }   (l = t,b <-> (-1)^l = +1,-1)
# =============================================================================

"""hbar^2/(2 m_e) in meV*Angstrom^2.  = 3.80998 eV*Angstrom^2"""
const HBAR2_2ME = 3809.98

"""2x2 rotation by `phi`"""
rot2(phi) = [cos(phi) -sin(phi); sin(phi) cos(phi)]

"""C3 of the paper: rotation by +2pi/3 (B10)."""
const C3M = rot2(2pi / 3)

C3op(v) = C3M * v

"""`n` points from `p` to `q` inclusive (endpoints shared)"""
function lineseg(p::AbstractVector, q::AbstractVector, n::Int)
    p = collect(float.(p))
    q = collect(float.(q))
    return [p .+ (q .- p) .* (j - 1) / (n - 1) for j in 1:n]
end

# -----------------------------------------------------------------------------
# MoireGeometry
# -----------------------------------------------------------------------------
struct MoireGeometry
    a0::Float64            # monolayer lattice constant
    theta::Float64         # twist angle (radians)
    aM::Float64            # moire lattice constant
    aM1::Vector{Float64}   # a_{M,1}
    aM2::Vector{Float64}   # a_{M,2}
    bM1::Vector{Float64}   # b_{M,1}
    bM2::Vector{Float64}   # b_{M,2}
    g::Vector{Vector{Float64}}    # g_1, g_2, g_3
    g2::Vector{Vector{Float64}}   # g_{2,1}, g_{2,2}, g_{2,3}
    q::Vector{Vector{Float64}}    # q_1, q_2, q_3
    q2::Vector{Vector{Float64}}   # q_{2,1}, q_{2,2}, q_{2,3}
    Gs::Vector{Vector{Float64}}   # all G_M in the hexagonal shell
    Gkey::Vector{NTuple{2,Int}}   # matching (m1,m2) integer labels
    Nshell::Int
    Rcut::Float64                 # radius of the +-K-valley plane-wave disk
    K_M::Vector{Float64}          # (b_1+b_2)/3        -- BZ corner at 30 deg
    Mp_M::Vector{Float64}         # midpoint(K_M, q_1)  -- BZ edge midpoint
end

"""
    moire_geometry(; a0 = 3.52, theta_deg = 3.89, Nshell = 5)

Build the moire geometry.  `Nshell` is the hexagonal truncation of the
reciprocal shell: |m1| <= N, |m2| <= N, |m1+m2| <= N  (same idea as
`moireband_PRL122.jl:79-101`, generalised).
"""
function moire_geometry(; a0::Real = 3.52, theta_deg::Real = 3.89, Nshell::Int = 5)
    theta = deg2rad(float(theta_deg))

    # (B1)
    aM = a0 / (2 * sin(theta / 2))

    # (B2)
    bM1 = 4pi / (sqrt(3) * aM) .* [1.0, 0.0]
    bM2 = 4pi / (sqrt(3) * aM) .* [0.5, sqrt(3) / 2]

    # (B1)
    aM1 = aM .* [sqrt(3) / 2, -0.5]
    aM2 = C3op(aM1)

    # g_i = C3^{i-1} b_{M,1}
    g = [copy(bM1)]
    push!(g, C3op(g[1]))
    push!(g, C3op(g[2]))

    # g_{2,1} = b_{M,1} + b_{M,2},  g_{2,i} = C3^{i-1} g_{2,1}
    g2 = [bM1 .+ bM2]
    push!(g2, C3op(g2[1]))
    push!(g2, C3op(g2[2]))

    # (B12) q_1 = 4pi/(3 a_0) * 2 sin(theta/2) * (0,1)
    q1 = 4pi / (3 * a0) * 2 * sin(theta / 2) .* [0.0, 1.0]
    q = [copy(q1)]
    push!(q, C3op(q[1]))
    push!(q, C3op(q[2]))

    # q_{2,1} = b_{M,1} + q_1
    q2 = [bM1 .+ q1]
    push!(q2, C3op(q2[1]))
    push!(q2, C3op(q2[2]))

    # hexagonal shell of G_M
    #
    # Empty typed vectors are spelled two ways in this file, on a rule: a nested
    # element type is written out (Vector{Vector{Float64}}()), a simple one uses
    # the sugar (NTuple{2,Int}[]).  The two are the same type -- T[] is exactly
    # Vector{T}() -- but the sugar for a *nested* element type drops a brace
    # level (Vector{Float64}[] is Vector{Vector{Float64}}, not
    # Vector{Vector{Float64}}()), so it is avoided for the G/Q lists.
    Gs = Vector{Vector{Float64}}()
    Gkey = NTuple{2,Int}[]
    for m1 in -Nshell:Nshell, m2 in -Nshell:Nshell
        if abs(m1 + m2) <= Nshell
            push!(Gs, m1 .* bM1 .+ m2 .* bM2)
            push!(Gkey, (m1, m2))
        end
    end

    # BZ corners: K_M = (b_1+b_2)/3 at 30 deg; the inequivalent corner is q_1
    # (at 90 deg).  The edge midpoint between them is the M point of the path.
    KM = (bM1 .+ bM2) ./ 3
    MpM = (KM .+ q1) ./ 2

    # Radius of the +-K-valley plane-wave set.
    #
    # The +-K-valley basis is {G_M +- q_1}.  It is *centred on +-q_1*, not on
    # G_M = 0, so truncating by the hex shell |m1|,|m2|,|m1+m2| <= N would NOT
    # be C3-invariant: C3 maps G + q_1 to C3 G + q_2 with q_2 = q_1 - b_{M,2},
    # i.e. to (C3 G - b_{M,2}) + q_1, and C3 G - b_{M,2} can fall outside the
    # shell.  We therefore truncate by |Q| instead, which is invariant under
    # both C3 (|C3 Q| = |Q|) and inversion (|-Q| = |Q|) by construction.  The
    # radius is inherited from the N-shell so that N still controls the size.
    Rcut = 0.0
    for G in Gs
        Rcut = max(Rcut, norm(G .+ q1), norm(G .- q1))
    end

    geo = MoireGeometry(a0, theta, aM, aM1, aM2, bM1, bM2,
                        g, g2, q, q2, Gs, Gkey, Nshell, Rcut, KM, MpM)
    check_geometry(geo)
    return geo
end

# -----------------------------------------------------------------------------
# self-checks of the geometry (run at construction time)
# -----------------------------------------------------------------------------
"""
    check_geometry(geo)

Assert every geometric identity the model relies on.  Throws on failure.
"""
function check_geometry(geo::MoireGeometry; tol::Float64 = 1e-10)
    aM1, aM2 = geo.aM1, geo.aM2
    bM1, bM2 = geo.bM1, geo.bM2

    # a_{M,i} . b_{M,j} = 2pi delta_ij
    @assert isapprox(dot(aM1, bM1), 2pi; atol = tol) "a1.b1 != 2pi"
    @assert isapprox(dot(aM2, bM2), 2pi; atol = tol) "a2.b2 != 2pi"
    @assert abs(dot(aM1, bM2)) < tol "a1.b2 != 0"
    @assert abs(dot(aM2, bM1)) < tol "a2.b1 != 0"

    # q_1 is a BZ corner: |q_1| == 4pi/(3 a_M) and q_1 is equidistant from 0
    # and b_{M,2} (i.e. it lies on the perpendicular bisector).
    nq1 = norm(geo.q[1])
    @assert isapprox(nq1, 4pi / (3 * geo.aM); atol = 1e-12) "|q1| wrong"
    @assert isapprox(nq1, norm(geo.q[1] .- bM2); atol = 1e-12) "q1 not on the bisector of b2"

    # the six +-q_i are distinct
    corners = Vector{Vector{Float64}}()
    for i in 1:3, s in (1.0, -1.0)
        push!(corners, s .* geo.q[i])
    end
    for i in 1:6, j in (i+1):6
        @assert norm(corners[i] .- corners[j]) > 1e-6 "degenerate BZ corners"
    end

    # K_M is a BZ corner of the *other* class: same length as q_1, and the angle
    # between them is 60 deg modulo the 120 deg point-group rotation (mod
    # 2pi/3 the difference must be +-pi/3, never 0 -- those are the two
    # inequivalent corner classes K and K').
    @assert isapprox(norm(geo.K_M), nq1; atol = 1e-12) "K_M is not a BZ corner"
    dang = atan(geo.K_M[2], geo.K_M[1]) - atan(geo.q[1][2], geo.q[1][1])
    @assert abs(abs(mod(dang, 2pi / 3)) - pi / 3) < 1e-10 "K_M and q_1 are equivalent corners"

    # 120-degree stars: sums vanish
    @assert norm(geo.g[1] .+ geo.g[2] .+ geo.g[3]) < tol "sum g_i != 0"
    @assert norm(geo.g2[1] .+ geo.g2[2] .+ geo.g2[3]) < tol "sum g_{2i} != 0"
    @assert norm(geo.q[1] .+ geo.q[2] .+ geo.q[3]) < tol "sum q_i != 0"
    @assert norm(geo.q2[1] .+ geo.q2[2] .+ geo.q2[3]) < tol "sum q_{2i} != 0"

    # C3 cyclicity
    for i in 1:3
        @assert norm(C3op(geo.g[i]) .- geo.g[mod1(i + 1, 3)]) < tol "C3 g_i cycle broken"
        @assert norm(C3op(geo.g2[i]) .- geo.g2[mod1(i + 1, 3)]) < tol "C3 g_{2i} cycle broken"
        @assert norm(C3op(geo.q[i]) .- geo.q[mod1(i + 1, 3)]) < tol "C3 q_i cycle broken"
        @assert norm(C3op(geo.q2[i]) .- geo.q2[mod1(i + 1, 3)]) < tol "C3 q_{2i} cycle broken"
    end

    # 3 q_1 is a moire reciprocal lattice vector (needed for w-term closure)
    @assert is_lattice_vec(geo, 3 .* geo.q[1]) "3 q_1 is not a G_M"
    # g_i, g_{2,i} are reciprocal lattice vectors
    for v in vcat(geo.g, geo.g2)
        @assert is_lattice_vec(geo, v) "g vector is not a G_M"
    end

    return true
end

"""is `v` an integer combination of b_{M,1}, b_{M,2}?"""
function is_lattice_vec(geo::MoireGeometry, v; tol::Float64 = 1e-9)
    n1, n2 = lattice_coeffs(geo, v)
    return abs(n1 - round(n1)) < tol && abs(n2 - round(n2)) < tol
end

"""coefficients (c1,c2) with v = c1 b_{M,1} + c2 b_{M,2}"""
function lattice_coeffs(geo::MoireGeometry, v)
    return (dot(v, geo.aM1) / (2pi), dot(v, geo.aM2) / (2pi))
end

# -----------------------------------------------------------------------------
# exact integer labelling of the Q sets
#
# Every Q is of the form  m1 b_{M,1} + m2 b_{M,2} + f q_1  with f = +-1, and
# 3 q_1 = -b_{M,1} + 2 b_{M,2}.  Hence  3Q = n1 b_{M,1} + n2 b_{M,2}  with
#     f = +1 :  (n1,n2) = (3m1 - 1, 3m2 + 2)
#     f = -1 :  (n1,n2) = (3m1 + 1, 3m2 - 2)
# so the family is read off exactly from the residues mod 3.  This makes every
# delta-function lookup exact integer arithmetic instead of floating point.
# -----------------------------------------------------------------------------

"""integer key (n1,n2) of Q = m1 b1 + m2 b2 + f q_1  <->  3Q = n1 b1 + n2 b2"""
qkey(m1::Int, m2::Int, f::Int) = (3m1 - f, 3m2 + 2f)

@inline key_add(k::NTuple{2,Int}, d::NTuple{2,Int}) = (k[1] + d[1], k[2] + d[2])
@inline vscale(d::NTuple{2,Int}, s::Int) = (s * d[1], s * d[2])

# --- K valley: the Q set {G_M + q_1} U {G_M - q_1} --------------------------
struct KValleyBasis
    Qs::Vector{Vector{Float64}}     # actual Q vectors
    f::Vector{Int}                  # +1 for {G_M + q_1}, -1 for {G_M - q_1}
    key::Vector{NTuple{2,Int}}      # exact integer key
    index::Dict{NTuple{2,Int},Int}  # key -> row/col index
end

@inline lookup(bas::KValleyBasis, k::NTuple{2,Int}) = get(bas.index, k, 0)

"""
    kvalley_basis(geo) -> KValleyBasis

The plane-wave set of (B23):  Q in {G_M + q_1} U {G_M - q_1}, truncated by
`|Q| <= geo.Rcut` (see the comment in `moire_geometry` on why the truncation is
by radius and not by shell index).  `f` records which family a state belongs
to; for eta = +1 the family f = +1 is the *top* layer (l = t, (-1)^l = +1) and
f = -1 the bottom layer (see (B26)).
"""
function kvalley_basis(geo::MoireGeometry)
    Qs = Vector{Vector{Float64}}()
    f = Int[]
    key = NTuple{2,Int}[]
    pad = geo.Nshell + 4                      # scan wide, then filter by radius
    for m1 in -pad:pad, m2 in -pad:pad, ff in (1, -1)
        Q = m1 .* geo.bM1 .+ m2 .* geo.bM2 .+ ff .* geo.q[1]
        norm(Q) <= geo.Rcut + 1e-9 || continue
        push!(Qs, Q)
        push!(f, ff)
        push!(key, qkey(m1, m2, ff))
    end
    index = Dict{NTuple{2,Int},Int}()
    for (i, k) in enumerate(key)
        @assert !haskey(index, k) "duplicate Q key -- the two families overlap"
        index[k] = i
    end

    # the truncation must be closed under both C3 and inversion, otherwise the
    # symmetry identities would fail for a purely artificial reason
    for (i, Q) in enumerate(Qs)
        for (name, QQ) in (("C3", C3op(Q)), ("C3^-1", rot2(-2pi / 3) * Q), ("-1", -Q))
            found = false
            for (j, Q2) in enumerate(Qs)
                if norm(QQ .- Q2) < 1e-8
                    found = true
                    break
                end
            end
            @assert found "K-valley basis is not closed under $(name): Q = $(Q)"
        end
    end
    return KValleyBasis(Qs, f, key, index)
end

# --- Gamma valley: triangular Q set {G_M} -----------------------------------
struct GammaValleyBasis
    Qs::Vector{Vector{Float64}}
    index::Dict{NTuple{2,Int},Int}
end

@inline lookup(bas::GammaValleyBasis, k::NTuple{2,Int}) = get(bas.index, k, 0)

function gammavalley_basis(geo::MoireGeometry)
    Qs = [copy(G) for G in geo.Gs]
    return GammaValleyBasis(Qs, Dict(k => i for (i, k) in enumerate(geo.Gkey)))
end

# -----------------------------------------------------------------------------
# high-symmetry paths through the moire BZ; the points are caller-supplied
# -----------------------------------------------------------------------------
struct KPath
    kpts::Vector{Vector{Float64}}
    nodes::Vector{Int}          # indices into kpts of the high-symmetry points, in order
end

"""
    bz_path(geo, points; nseg = fill(40, length(points) - 1)) -> KPath

Piecewise-linear k-path visiting `points` in order.  `points` is supplied by the
caller -- nothing about the path is fixed here, so an arbitrary high-symmetry
route can be requested.  For the route of the paper's Fig. 4 write

    bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])

`nseg[j]` is the number of k-points on segment `j` (the straight line
`points[j] -> points[j+1]`), endpoints included; a shared endpoint is stored only
once, so `kpts` has `sum(nseg) - (length(points) - 1)` entries.  `nodes[j]` is the
index in `kpts` of `points[j]`, i.e. `kpts[nodes[j]] == points[j]`.
"""
function bz_path(geo::MoireGeometry, points::AbstractVector{<:AbstractVector};
                 nseg::AbstractVector{<:Integer} = fill(40, length(points) - 1))
    @assert length(points) >= 2 "need at least two high-symmetry points"
    @assert length(nseg) == length(points) - 1 "nseg must have one entry per segment"
    @assert all(>(1), nseg) "each segment needs at least two points"

    kpts = Vector{Vector{Float64}}()
    nodes = Int[]
    for j in 1:(length(points) - 1)
        seg = lineseg(points[j], points[j + 1], nseg[j])
        for (t, p) in enumerate(seg)
            t == 1 && !isempty(kpts) && continue  # node already stored by the previous segment
            push!(kpts, p)
        end
        push!(nodes, length(kpts))                # the segment's last point is points[j+1]
    end
    pushfirst!(nodes, 1)                          # points[1] is kpts[1]
    return KPath(kpts, nodes)
end
