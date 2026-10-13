# =============================================================================
# 03_selftest.jl -- validate the whole chain before drawing anything
#
# Run in the *default* environment (so the already-installed Plots is visible):
#
#     julia scripts/03_selftest.jl
#
# The judgement criteria are the paper's own symmetry identities, not a visual
# comparison with Fig. 4.  Everything asserted here is a matrix identity between
# two independently constructed operators, an exact degeneracy the paper
# predicts, or an agreement between the numerics and an independent analytic
# formula from the appendices.
# =============================================================================

using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD

const KGEN = [0.0173, -0.0091]     # a generic k, off every high-symmetry line

println()
println("Julia $(VERSION)   ", Sys.MACHINE)
println("θ = $(THETA_DEFAULT_DEG)°, a0 = $(A0_MOTE2) Å, Nshell = 5")
println()

geo = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = 5)
@printf("a_M = %.4f Å, |b_M1| = %.6f Å⁻¹, |q1| = %.6f Å⁻¹, |K_M| = %.6f Å⁻¹\n",
        geo.aM, norm(geo.bM1), norm(geo.q[1]), norm(geo.K_M))
let path = bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])
    @printf("path nodes: Γ_M@%d  K_M@%d  M_M@%d  Γ_M@%d  (nk = %d)\n",
            path.nodes..., length(path.kpts))
end
println()

# --- Part A: symmetry / consistency identities -------------------------------
results = symmetry_report(geo; eps = 0.0)
ok = print_report(results)

# --- Part B: the two sources of h_QQ'(k) must agree element by element --------
# "Main source": the analytic Fourier components of the *real-space* potentials
# (B21)/(B38)/(C4)/(C13), assembled in kmodels.jl.  "Cross source": the printed
# *momentum-space* equations (B27)/(C8)/(C18), transcribed verbatim in
# transcribe.jl.  The two were written by different routes -- integer-key
# bookkeeping over a Fourier table vs. a float delta-function scan over every
# pair of basis states -- so exact agreement validates both.  (App. B 3 prints
# no momentum-space Gamma form for AA stacking, so that model has no second
# source; its symmetry identities are checked in Part A instead.)
println()
println("Cross-source: kmodels.jl (real-space Fourier components) vs transcribe.jl (printed δ-form)")
resid = Float64[]
for (model, prm) in ((:AA_K, AA_K_FH), (:AA_K, AA_K_FHSH), (:AB_K, AB_K),
                     (:AB_Gamma, AB_GAMMA))
    bas = build_basis(geo, model)
    for eta in (model in (:AA_K, :AB_K) ? (1, -1) : (1,))
        Hk = hamiltonian(geo, bas, prm, model, KGEN; eta = eta, convention = :paper)
        Ht = transcribe_hamiltonian(geo, bas, prm, model, KGEN; eta = eta)
        Hr = hamiltonian(geo, bas, prm, model, KGEN; eta = eta, convention = :transposed)
        v1 = norm(Hk - Ht)                       # the two sources, same reading
        v2 = norm(Ht - transpose(Hr))            # and the psi -> -psi reading
        @printf("  %-10s %-22s η=%+d   ‖kmodels − transcribe‖ = %9.2e   ‖transcribe − transp H‖ = %9.2e\n",
                model, label(prm), eta, v1, v2)
        push!(resid, v1, v2)
    end
end
@printf("  worst residual: %.2e\n", maximum(resid))
ok &= maximum(resid) < 1e-10

# --- Part C: ε ≠ 0 must break the effective inversion / mirror-T symmetries ---
# (B19)/(B40)/(C15): the displacement field is (-1)^l ε/2, which flips sign
# under a layer exchange.  Every identity whose operator swaps the layers must
# therefore fail at ε ≠ 0 -- and that is the *only* lever that exposes the layer
# structure of the Γ-valley permutations, because at ε = 0 the Γ models are
# layer-exchange symmetric on their own and pass either way.  (The ε term is
# diagonal in the layer, so it cannot move the Γ-valley flat bands either: that
# is the paper's remark at the end of App. C 3.)
println()
println("ε = 10 meV: everything that flips the layer must FAIL (B19)/(B40)/(C15)")
res_eps = symmetry_report(geo; eps = 10.0)
invs = filter(c -> occursin("effective inversion", c.name) ||
                   occursin("effective K", c.name) ||
                   occursin("C2yT", c.name) || occursin("C2xT", c.name), res_eps)
for c in invs
    okstr = occursin("AB±K", c.name) && occursin("effective inversion", c.name) ?
            "(never a symmetry of AB)" :
            (c.passed ? "UNEXPECTED PASS — layer structure not probed" :
                        "ok (broken as predicted)")
    @printf("  %-46s %10.2e  %s\n", c.name, c.value, okstr)
end

# --- Part D: shell convergence of the ±K low-energy bands --------------------
# The truncation must be converged before any number is quoted: the plane-wave
# set grows with Nshell, so the low bands can only move *down* (variational).
# A zero here would mean the extra states are simply not being used, so the
# basis dimension is printed alongside the energy change.
println()
println("Shell convergence: top 6 valence bands vs Nshell, at Γ_M")
for (model, prm, n) in ((:AA_K, AA_K_FH, 6), (:AA_K, AA_K_FHSH, 6), (:AB_K, AB_K, 6))
    Es = Vector{Vector{Float64}}()
    dims = Int[]
    for N in (4, 5, 6)
        g = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = N)
        b = build_basis(g, model)
        E = diagonalize_point(g, b, prm, model, [0.0, 0.0])
        push!(Es, E[(end - n + 1):end])
        push!(dims, b isa GammaBasis ? b.n : length(b.Qs))
    end
    @printf("  %-10s %-22s  dim %3d/%3d/%3d   N=4→5: %9.3e   N=5→6: %9.3e meV\n",
            model, label(prm), dims..., maximum(abs.(Es[1] .- Es[2])),
            maximum(abs.(Es[2] .- Es[3])))
end

# --- Part E: the Γ-valley flat bands are harmonic-oscillator levels ----------
# App. C 3 argues that -V_Γ(r) confines the top Γ-valley states to |r| ~ 3.7 Å
# inside a 52 Å moiré cell, so they are the (n_x,n_y) eigenstates of an array of
# decoupled harmonic oscillators.  (C24)/(C25) give the harmonic spectrum: ħω =
# 56.76 meV, the (0,0) level single (2-fold with spin) with spinless C3 = 1, and
# the (1,0),(0,1) doublet (4-fold with spin) with C3 = e^{±i2π/3}.  The exact
# lattice model must reproduce that structure, with a splitting that differs
# from ħω by the quartic anharmonicity -- an independent analytic prediction
# (parameters.jl:anharmonic_shift) that is checked numerically here.
println()
println("Γ-valley flat bands (AB): exact lattice vs (C24) harmonic oscillator")
ok_gamma = let g7 = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = 7)
    bG = build_basis(g7, :AB_Gamma)
    Hs = ham_AB_Gamma(g7, bG, AB_GAMMA, [0.0, 0.0])         # spinless, 2 layers
    Pc3 = Matrix{Float64}(I, size(Hs, 1), size(Hs, 1))[c3_perm_spinless(bG), :]
    lev = top_multiplets(Hs, Pc3; nlev = 2)
    l0, l1 = lev[1], lev[2]
    hw = hbar_omega(g7, AB_GAMMA)
    sh = anharmonic_shift(g7, AB_GAMMA)
    resid_g = abs((l0.E - l1.E) - (hw + sh))
    @printf("  (0,0):        %d state,   spinless C3 = %.6f%+.6fi\n",
            l0.mult, real(l0.c3[1]), imag(l0.c3[1]))
    @printf("  (1,0),(0,1):  %d states,  spinless C3 = %s\n", l1.mult,
            join([@sprintf("%.4f%+.4fi", real(x), imag(x)) for x in l1.c3], "  "))
    @printf("  splitting: exact lattice %.3f meV   ħω (C25) %.3f   ħω + anharmonic %.3f\n",
            l0.E - l1.E, hw, hw + sh)
    @printf("  |exact − (ħω + first-order quartic)| = %.4f meV\n", resid_g)
    good = l0.mult == 1 && abs(l0.c3[1] - 1) < 1e-8 && l1.mult == 2 && resid_g < 0.25

    # The splitting is slowly convergent in Nshell (the Γ-valley oscillator is
    # far more localized than the ±K states), so show the sequence: it is a
    # single-Γ_M diagonalisation, hence cheap.
    @printf("  splitting vs Nshell:")
    seq = Float64[]
    for N in (4, 5, 6, 7, 8, 9)
        g = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = N)
        b = build_basis(g, :AB_Gamma)
        lv = top_multiplets(ham_AB_Gamma(g, b, AB_GAMMA, [0.0, 0.0]),
                            Matrix{Float64}(I, 2length(b.G.Qs), 2length(b.G.Qs))[
                                c3_perm_spinless(b), :]; nlev = 2)
        push!(seq, lv[1].E - lv[2].E)
        @printf(" %s=%.2f", N, seq[end])
    end
    @printf("  meV   (N=8→9: %.4f)\n", abs(seq[5] - seq[6]))
    good &= abs(seq[5] - seq[6]) < 0.05

    # "Extremely flat": localized on a 3.7 Å scale, so the dispersion over the
    # 52 Å moiré cell is a small fraction of ħω.
    path4 = bz_path(g7, [[0.0, 0.0], g7.K_M, g7.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])
    b0 = [diagonalize_point(g7, bG, AB_GAMMA, :AB_Gamma, kk)[(end-1):end]
          for kk in path4.kpts]
    b1 = [diagonalize_point(g7, bG, AB_GAMMA, :AB_Gamma, kk)[(end-3):(end-2)]
          for kk in path4.kpts]
    bw(v) = maximum(v) - minimum(v)
    w0 = max(bw([x[1] for x in b0]), bw([x[2] for x in b0]))
    w1 = max(bw([x[1] for x in b1]), bw([x[2] for x in b1]))
    @printf("  bandwidth over the MBZ: (0,0) %.4f meV,  (1,0),(0,1) %.4f meV\n", w0, w1)
    good &= max(w0, w1) < 0.05

    # E_Γ is not tabulated anywhere; the default is *calibrated* by the paper's
    # "30 meV below the ±K valence band maximum" statement (Sec. III.C).  The
    # calibration must be reproducible.
    rec = E0_relative_to_vbm(g7, AB_GAMMA, AB_K; target = -30.0)
    @printf("  E_Γ: default %10.4f meV, recomputed from the VBM rule %10.4f meV (Δ = %.1e)\n",
            payload(AB_GAMMA).E0, payload(rec).E0,
            abs(payload(rec).E0 - payload(AB_GAMMA).E0))
    good &= abs(payload(rec).E0 - payload(AB_GAMMA).E0) < 0.05
    good
end
ok &= ok_gamma

# --- Part F: AA ±K band structure, sanity numbers ----------------------------
println()
println("AA ±K (FH+SH), top valence bands at the high-symmetry points (meV)")
geo5 = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = 5)
bas = build_basis(geo5, :AA_K)
for (name, kk) in (("Γ_M", [0.0, 0.0]), ("K_M", geo5.K_M), ("M_M", geo5.Mp_M))
    E = diagonalize_point(geo5, bas, AA_K_FHSH, :AA_K, kk)
    @printf("  %-4s E4..E1 = %8.3f %8.3f %8.3f %8.3f\n", name,
            E[end-3], E[end-2], E[end-1], E[end])
end

# =============================================================================
# Part G: the parameter bookkeeping (ParamSet / ParamRef / provenance guard)
#
# All of Parts A-F above are the *physics*; this part is the *record keeping*
# around it.  It has two jobs.
#
#   (1) Fence: every number that used to be a bare literal is still that exact
#       number, and the whole band structure is bit-for-bit unchanged.  The
#       provenance metadata must be physically inert -- if this section moves a
#       checksum, the refactor changed the model, not just the paperwork.
#
#   (2) Guard: Table I / Table II are fits made at ONE material and ONE twist
#       angle (PRB 109, 205121 (2024), theta = 3.89 deg, MoTe2).  Using them
#       elsewhere produces smooth, plausible, wrong bands, so `check_parameters`
#       must throw -- and each of the two escape hatches must work, because
#       crossing a boundary on purpose is legitimate.
# =============================================================================
println()
println("Part G: parameter records, provenance guard, and the refactor's inertness")
ok_g = let
    g = true

# --- G1. the fitted numbers are bit-identical to the pre-ParamSet literals ----
frozen = [
    (AA_K_FH,      KParams,   (0.60, 16.5, deg2rad(-105.9), -18.8,  0.0,   0.0)),
    (AA_K_FHSH,    KParams,   (0.62,  7.94, deg2rad(-88.43), -10.77, 20.00, 10.21)),
    (AB_K,         KParams,   (0.62, 53.0, deg2rad(-56.0),   0.0,   0.0,   0.0)),
    (AB_GAMMA,     GammaParams, (10.0, 72.0, 0.0, 300.0, -640.1296)),
    (AA_GAMMA,     GammaParams, (10.0, 72.0, 0.0, 300.0, -817.3836)),
]
for (ps, T, want) in frozen
    got = T === KParams ?
          (payload(ps).mstar, payload(ps).V, payload(ps).psi, payload(ps).w,
           payload(ps).V2, payload(ps).w2) :
          (payload(ps).mstar, payload(ps).V, payload(ps).psi, payload(ps).w,
           payload(ps).E0)
    same = all(map(==, got, want))
    @printf("  %-12s fields == frozen literals    %s\n", label(ps), same ? "yes" : "NO: $got")
    g &= same
end

# ...and the geometry they claim must stay the paper's geometry.  The five
# records express `a0`/`theta_deg` through the package defaults, so editing a
# default (`geometry.jl`) would silently re-point every one of them -- and the
# guard, whose whole job is to catch a parameters/geometry mismatch, would go
# quiet instead of loud.  Pinning the defaults here forces a deliberate choice.
geom_ok = A0_DEFAULT == 3.52 && THETA_DEFAULT_DEG == 3.89
@printf("  package defaults still the paper's geometry (a0 = 3.52, theta = 3.89)  %s\n",
        geom_ok ? "yes" : "NO -- the five records below now claim the new geometry")
g &= geom_ok

# --- G2. registry: derived from the constants, not a second source of truth ---
reg_ok = length(PAPER_FITS) == 5
for (key, ps) in PAPER_FITS
    r = paramref(ps)
    reg_ok &= (r.stacking, r.valley, r.harmonics) == key
    reg_ok &= paper_fit(key...) === ps                 # identity, not a copy
end
@printf("  registry: %d keys, each ref matches its key, paper_fit returns the same object  %s\n",
        length(PAPER_FITS), reg_ok ? "yes" : "NO")
throws(f) = try (f(); false) catch; true end
nofit = throws(() -> paper_fit(:AB, :K, :FHSH))        # no such fit in the paper
@printf("  paper_fit on an absent key throws  %s\n", nofit ? "yes" : "NO")
g &= reg_ok && nofit

# --- G3. the guard fires exactly where it must -------------------------------
# Every recorded set is valid at its own fit geometry...
fit_geo = moire_geometry(a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG)
passes = true
for (model, ps) in ((:AA_K, AA_K_FH), (:AA_K, AA_K_FHSH), (:AB_K, AB_K),
                    (:AA_Gamma, AA_GAMMA), (:AB_Gamma, AB_GAMMA))
    passes &= check_parameters(fit_geo, ps; model = model) === true
end
@printf("  all five sets accepted at their own fit geometry (a0 = %g, theta = %g deg)  %s\n",
        A0_MOTE2, THETA_DEFAULT_DEG, passes ? "yes" : "NO")

# ...and nowhere else.  These three are the failure modes that matter.
other_mat = moire_geometry(a0 = 3.282, theta_deg = THETA_DEFAULT_DEG)   # a WSe2-like a0
other_ang = moire_geometry(a0 = A0_MOTE2, theta_deg = 5.0)
cases = [
    ("material mismatch   a0 = 3.282", () -> check_parameters(other_mat, AB_GAMMA; model = :AB_Gamma)),
    ("angle mismatch      theta = 5 deg", () -> check_parameters(other_ang, AA_K_FHSH; model = :AA_K)),
    ("stacking/valley mismatch (AB_K at :AA_K)", () -> check_parameters(fit_geo, AB_K; model = :AA_K)),
]
for (what, f) in cases
    fired = throws(f)
    @printf("  %-42s throws  %s\n", what, fired ? "yes" : "NO")
    g &= fired
end

# --- G4. the two escape hatches, and that they are reversible ---------------
esc1 = check_parameters(other_ang, exploratory(AA_K_FHSH); model = :AA_K) === true
enforce_param_checks!(false)
esc2 = check_parameters(other_ang, AA_K_FHSH; model = :AA_K) === true
enforce_param_checks!(true)
restored = throws(() -> check_parameters(other_ang, AA_K_FHSH; model = :AA_K))
@printf("  escape hatches: exploratory(ps) %s, enforce_param_checks!(false) %s, restored %s\n",
        esc1 ? "ok" : "NO", esc2 ? "ok" : "NO", restored ? "ok" : "NO")
g &= esc1 && esc2 && restored

# the provenance survives the escape hatch: same payload, claim withdrawn
esc_ok = payload(exploratory(AA_K_FHSH)) === payload(AA_K_FHSH) &&
         paramref(exploratory(AA_K_FHSH)).provenance === EXPLORATORY &&
         paramref(AA_K_FHSH).provenance === FITTED
@printf("  exploratory() keeps the payload and withdraws the geometry claim  %s\n",
        esc_ok ? "yes" : "NO")
g &= esc_ok

# --- G5. Material == the raw a0, and the calibrated set re-stamps itself -----
geo_mat = moire_geometry(MOTE2; theta_deg = THETA_DEFAULT_DEG)
geo_raw = moire_geometry(a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG)
mat_ok = geo_mat.a0 == geo_raw.a0 && geo_mat.aM == geo_raw.aM &&
         geo_mat.bM1 == geo_raw.bM1 && geo_mat.q == geo_raw.q &&
         geo_mat.K_M == geo_raw.K_M
@printf("  moire_geometry(MOTE2) == moire_geometry(a0 = 3.52)  %s\n", mat_ok ? "yes" : "NO")
g &= mat_ok

# The calibration is a *number*, so it is only reproducible at the truncation
# the stored value was computed with: the Gamma-valley flat band needs the
# Nshell = 7 basis (Part E), and 0.37 meV of its energy is still missing at
# Nshell = 5.  Note this is not a geometry mismatch -- a0 and theta, the two
# things the guard watches, do not depend on the cutoff.
geo7 = moire_geometry(MOTE2; theta_deg = THETA_DEFAULT_DEG, Nshell = 7)
rec = E0_relative_to_vbm(geo7, AB_GAMMA, AB_K; target = -30.0)
rec5 = E0_relative_to_vbm(geo_mat, AB_GAMMA, AB_K; target = -30.0)
rec_ref = paramref(rec)
cal_ok = rec_ref.provenance === CALIBRATED &&
         rec_ref.calibrated_at == (geo7.a0, rad2deg(geo7.theta)) &&
         abs(payload(rec).E0 - (-640.1296)) < 0.05
@printf("  E0_relative_to_vbm: provenance %s, calibrated_at == (a0, theta) of the geometry %s, E_Γ = %.4f meV (Nshell = 5: %.4f)  %s\n",
        rec_ref.provenance, rec_ref.calibrated_at == (geo7.a0, rad2deg(geo7.theta)) ? "yes" : "NO",
        payload(rec).E0, payload(rec5).E0, cal_ok ? "ok" : "NO")
g &= cal_ok
# and the re-stamped record is refused at any other angle
cal_guard = throws(() -> check_parameters(other_ang, rec; model = :AB_Gamma))
@printf("  ... and that calibrated record is refused at theta = 5 deg  %s\n",
        cal_guard ? "yes" : "NO")
g &= cal_guard

# --- G6. the reference fence: the bands themselves did not move --------------
# `ref` is metadata; it must not touch a single floating-point operation.
baseline = [(:AB_K, AB_K, 4.26562296284693778e+07),
            (:AA_K, AA_K_FH, 4.40515890905550867e+07),
            (:AA_K, AA_K_FHSH, 4.26288541207873896e+07)]
for (model, ps, want) in baseline
    got = sum(abs, diagonalize_path(fit_geo, model, ps).energies)
    same = isapprox(got, want; rtol = 1e-12)
    @printf("  Σ|E| along ΓKMΓ  %-10s %-18s %.17e  %s\n", model, label(ps), got,
            same ? "unchanged" : "MOVED (frozen: $(want))")
    g &= same
end
vbm_now = vbm(fit_geo, :AB_K, AB_K)
vbm_ok = isapprox(vbm_now, 6.60526527534205883e+01; rtol = 1e-12)
@printf("  vbm(:AB_K) %.17e  %s\n", vbm_now, vbm_ok ? "unchanged" : "MOVED")
g &= vbm_ok
g &= isapprox(hbar_omega(fit_geo, AB_GAMMA), 5.67614326631262074e+01; rtol = 1e-12)
g &= isapprox(anharmonic_shift(fit_geo, AB_GAMMA), -1.86450245252928859e+00; rtol = 1e-12)

    g
end
ok &= ok_g

println()
if ok
    println("SELF-TEST PASSED")
else
    println("SELF-TEST FAILED")
    exit(1)
end
