# =============================================================================
# bands.jl -- diagonalise h_QQ'(k), and the E_Gamma calibration rule
#
# There is no band-structure object here.  A k-path is whatever list of k
# points the caller assembles with `bz_path`, and `diagonalize_point` is mapped
# over it; the eigenvalues are then the caller's to slice, shift and plot.
# Nothing in this file knows about the paper's Fig. 4 route -- the only two
# things it needs from a k-set are a maximum (for `vbm`) and a rule for fixing
# the free offset E_Gamma.
# =============================================================================

"""eigenvalues of h(k), ascending"""
function diagonalize_point(geo::MoireGeometry, bas, p, model::Symbol, k;
                           eta::Int = 1, eps::Real = 0.0, convention::Symbol = :paper)
    H = hamiltonian(geo, bas, p, model, k; eta = eta, eps = eps, convention = convention)
    return sort!(real.(eigvals(Hermitian(H))))
end

"""
    vbm(geo, bas, p, model, kpts; eta = 1, eps = 0.0) -> Float64
    vbm(geo, bas, p, model, path::KPath; kw...)

The valence band maximum of `model`, i.e. the highest eigenvalue anywhere in
the k-set `kpts` (a vector of k points, or a `KPath` from `bz_path`).

`kpts` is a *sampling*, not the Brillouin zone: this is the maximum over the
points given, so a coarse set can miss the true VBM.  Which set is right is a
question about the model, not about this function, which is why there is no
default -- assemble the route you mean with `bz_path` and pass it.
"""
function vbm(geo::MoireGeometry, bas, p, model::Symbol,
             kpts::AbstractVector{<:AbstractVector}; kw...)
    isempty(kpts) && throw(ArgumentError("vbm: the k-set is empty"))
    return maximum(diagonalize_point(geo, bas, p, model, k; kw...)[end] for k in kpts)
end

vbm(geo::MoireGeometry, bas, p, model::Symbol, path::KPath; kw...) =
    vbm(geo, bas, p, model, path.kpts; kw...)

"""
    valence_bands(geo, bas, p, model, path; nbands, eta = 1, eps = 0.0)
        -> (d, nodes, E)

The top `nbands` valence bands of `model` along `path`, ready to plot: `d` is
the cumulative path length (`path_lengths`), `nodes` the indices of the
high-symmetry points, and `E` an `nk × nbands` matrix of eigenvalues in meV,
ascending in each row.  Nothing is shifted -- the caller picks the energy zero
(with `vbm`, say) and subtracts it.

A path is a sampling, so there is no default: assemble it with `bz_path`.
"""
function valence_bands(geo::MoireGeometry, bas, p, model::Symbol, path::KPath;
                       nbands::Int, eta::Int = 1, eps::Real = 0.0)
    E = reduce(hcat, [diagonalize_point(geo, bas, p, model, k; eta = eta, eps = eps)
                      for k in path.kpts])'
    n = size(E, 2)
    nbands <= n || throw(ArgumentError("valence_bands: asked for $nbands bands, " *
                                       "the basis has $n"))
    return (d = path_lengths(path.kpts), nodes = path.nodes,
            E = Matrix(E[:, (n - nbands + 1):n]))
end

"""
    E0_relative_to_vbm(geo, p_gamma, p_K, kpts; model_gamma = :AB_Gamma,
                       model_K = :AB_K, target = -30.0, eta = 1, eps = 0.0)
        -> ParamSet

`E_Gamma` of (C19) is not tabulated anywhere in the paper -- it is the rigid
offset between the inactive Gamma valley and the active +-K valleys (App. B 3),
fixed in practice by requiring the first ultra-flat Gamma band to lie `target`
meV below the +-K valence band maximum.  Returns `p_gamma` with the `E0` that
does exactly that, measured on the k-set `kpts` (so the answer is a property of
that sampling; the paper reads it off the Fig. 4(c) route, which the caller
supplies with `bz_path`).

The result is only meaningful at the geometry it was calibrated at, so its
provenance record is re-stamped with `geo`'s (a0, theta) and `CALIBRATED`, and
`check_parameters` will stop it from being carried to another geometry.  Both
`p_gamma` and `p_K` must themselves be valid at `geo` -- the paper's V, w and m*
are single-angle fits, so calibrating at a new angle first needs a refit (or an
explicit `exploratory`).
"""
function E0_relative_to_vbm(geo::MoireGeometry, p_gamma::ParamSet, p_K,
                            kpts::AbstractVector{<:AbstractVector};
                            model_gamma::Symbol = :AB_Gamma, model_K::Symbol = :AB_K,
                            target::Real = -30.0, eta::Int = 1, eps::Real = 0.0)
    gm = payload(p_gamma)
    # the Gamma Hamiltonian is linear in E0, so one probe at E0 = 0 fixes it
    probe = ParamSet(p_gamma.ref, GammaParams(gm.mstar, gm.V, gm.psi, gm.w, 0.0, gm.label))
    E0 = vbm(geo, build_basis(geo, model_K), p_K, model_K, kpts;
             eta = eta, eps = eps) + target -
         vbm(geo, build_basis(geo, model_gamma), probe, model_gamma, kpts;
             eta = eta, eps = eps)
    out = GammaParams(gm.mstar, gm.V, gm.psi, gm.w, E0, gm.label)
    ref = paramref(p_gamma)
    a0, theta_deg = geo.a0, rad2deg(geo.theta)
    return ParamSet(ParamRef(material = ref.material, a0 = a0, theta_deg = theta_deg,
                             stacking = ref.stacking, valley = ref.valley,
                             harmonics = ref.harmonics, provenance = CALIBRATED,
                             source = ref.source * "; E_Gamma calibrated here by the " *
                                      "'$(target) meV below the +-K VBM' rule",
                             calibrated_at = (a0, theta_deg)), out)
end

E0_relative_to_vbm(geo::MoireGeometry, p_gamma::ParamSet, p_K, path::KPath; kw...) =
    E0_relative_to_vbm(geo, p_gamma, p_K, path.kpts; kw...)
