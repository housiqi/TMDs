# =============================================================================
# bands.jl -- diagonalise h_QQ'(k) along a k-path
# =============================================================================

"""energies, in meV (ascending), plus the number of plane-wave states used"""
struct BandData
    energies::Matrix{Float64}     # nk x nstates, ascending in energy
    kdist::Vector{Float64}        # cumulative |dk| along the path
    nodes::Vector{Int}            # indices of Gamma, K, M, Gamma
    labels::Vector{String}
    nstates::Int
    model::Symbol
    eta::Int
end

"""eigenvalues of h(k), ascending"""
function diagonalize_point(geo::MoireGeometry, bas, p, model::Symbol, k;
                           eta::Int = 1, eps::Real = 0.0, convention::Symbol = :paper)
    H = hamiltonian(geo, bas, p, model, k; eta = eta, eps = eps, convention = convention)
    return sort!(real.(eigvals(Hermitian(H))))
end

"""
    diagonalize_path(geo, model, p; path = bz_path(geo), eta = 1, eps = 0.0,
                     convention = :paper) -> BandData
"""
function diagonalize_path(geo::MoireGeometry, model::Symbol, p;
                          path::KPath = bz_path(geo),
                          eta::Int = 1, eps::Real = 0.0,
                          convention::Symbol = :paper)
    bas = build_basis(geo, model)
    nk = length(path.kpts)
    E = reduce(hcat, [diagonalize_point(geo, bas, p, model, k;
                                        eta = eta, eps = eps, convention = convention)
                      for k in path.kpts])'      # nk x nstates
    # cumulative distance along the path
    d = zeros(nk)
    for i in 2:nk
        d[i] = d[i-1] + norm(path.kpts[i] .- path.kpts[i-1])
    end
    return BandData(E, d, path.nodes, ["Γ_M", "K_M", "M_M", "Γ_M"],
                    size(E, 2), model, eta)
end

"""the top `n` bands (highest valence states) of a BandData"""
top_bands(bd::BandData, n::Int) = bd.energies[:, (bd.nstates - n + 1):bd.nstates]

"""
    rigid_shift!(bd, ref)

Shift every band down by the valence band maximum so that the top band touches
0 meV, as in Fig. 4 of the paper (the paper plots energies relative to the VBM;
the absolute offset is a gauge choice of the two-center model).
"""
function rigid_shift!(bd::BandData)
    bd.energies .-= maximum(bd.energies)
    return bd
end

"""shift by an explicit energy instead"""
function shifted(bd::BandData, e0::Real)
    BandData(bd.energies .- e0, bd.kdist, bd.nodes, bd.labels,
             bd.nstates, bd.model, bd.eta)
end

"""the valence band maximum of `model` along the path (highest eigenvalue anywhere)"""
function vbm(geo::MoireGeometry, model::Symbol, p; path::KPath = bz_path(geo),
             eta::Int = 1, eps::Real = 0.0)
    bas = build_basis(geo, model)
    return maximum(diagonalize_point(geo, bas, p, model, k;
                                     eta = eta, eps = eps)[end] for k in path.kpts)
end

"""copy of a `GammaParams` with a different `E_Gamma`"""
with_E0(p::GammaParams, E0::Real) = GammaParams(p.mstar, p.V, p.psi, p.w,
                                                float(E0), p.label)

"""
    E0_relative_to_vbm(geo, p_gamma, p_K; model_gamma = :AB_Gamma,
                       model_K = :AB_K, target = -30.0, path = bz_path(geo))
        -> GammaParams

`E_Gamma` of (C19) is not tabulated anywhere in the paper -- it is the rigid
offset between the inactive Gamma valley and the active +-K valleys (App. B 3),
fixed in practice by requiring the first ultra-flat Gamma band to lie `target`
meV below the +-K valence band maximum.  Returns `p_gamma` with the `E0` that
does exactly that, at whatever twist angle `geo` describes.  The default
`target = -30` is the paper's Fig. 4(c) statement for AB stacking.
"""
function E0_relative_to_vbm(geo::MoireGeometry, p_gamma::GammaParams, p_K;
                            model_gamma::Symbol = :AB_Gamma, model_K::Symbol = :AB_K,
                            target::Real = -30.0, path::KPath = bz_path(geo))
    # the Gamma Hamiltonian is linear in E0, so one probe at E0 = 0 fixes it
    probe = with_E0(p_gamma, 0.0)
    return with_E0(p_gamma, vbm(geo, model_K, p_K; path = path) + target -
                           vbm(geo, model_gamma, probe; path = path))
end
