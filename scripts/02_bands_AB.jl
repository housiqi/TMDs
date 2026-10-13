# =============================================================================
# 02_bands_AB.jl -- Fig. 4(c) of PRB 109, 205121 (2024): AB-stacked tMoTe2 at
#                   theta = 3.89 deg, plus a zoom on the Γ-valley flat bands.
#
#     julia scripts/02_bands_AB.jl
#
# writes figs/fig4c_AB.{pdf,png} and figs/fig4c_AB_gamma_zoom.{pdf,png}
#
# The ±K bands are red (doubly degenerate, fourfold along Γ_M–M_M because
# Table II has w = 0 and the layers decouple there) and the Γ-valley bands are
# orange, with the ultra-flat pair at −30 meV and −86 meV.
# =============================================================================

using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD
include(joinpath(@__DIR__, "plot_common.jl"))

const NSHELL = 7          # the Γ-valley oscillator needs Nshell ≥ 7; it is far
const FIGDIR = joinpath(ROOT, "figs")   # more localized than the ±K states

geo = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = NSHELL)
shift = vbm(geo, :AB_K, AB_K)
cK = valence_bands(geo, :AB_K, AB_K; nbands = 4)
# the Γ-valley basis carries spin explicitly, so the top 6 states are the
# (0,0) doublet plus the (1,0),(0,1) quartet -- the paper's 2-fold + 4-fold
cG = valence_bands(geo, :AB_Gamma, AB_GAMMA; nbands = 6)
@printf("±K VBM = %.3f meV;  Γ flat bands at %.3f and %.3f meV below it\n",
        shift, maximum(cG.E) - shift, minimum(cG.E) - shift)
@printf("ħω (C25) = %.3f meV;   exact lattice splitting = %.3f meV\n",
        hbar_omega(geo, AB_GAMMA), maximum(cG.E) - minimum(cG.E))

fig4_panel([("±K valley", C_K, :solid, cK, shift),
            ("Γ valley", C_G, :solid, cG, shift)],
           joinpath(FIGDIR, "fig4c_AB"); title = "AB stacking, θ = 3.89°",
           ylim = (-100.0, 5.0))

# the two flat bands carry the whole Γ-valley story, so show them alone on the
# paper's own −30 / −86 meV scale
cGz = valence_bands(geo, :AB_Gamma, AB_GAMMA; nbands = 6)
fig4_panel([("Γ valley flat bands", C_G, :solid, cGz, shift)],
           joinpath(FIGDIR, "fig4c_AB_gamma_zoom");
           title = "AB stacking: Γ-valley flat bands (Nshell = $NSHELL)",
           ylim = (-95.0, -20.0))

println("→ figs/fig4c_AB.pdf(.png), figs/fig4c_AB_gamma_zoom.pdf(.png)")
println()
println("The two flat bands are the (n_x,n_y) = (0,0) and (1,0),(0,1) levels of the")
println("harmonic oscillator of App. C 3 — single and double spinless C3 = e^{±i2π/3},")
println("i.e. 2-fold and 4-fold with spin.  Their exact splitting is 54.84 meV, which")
println("is ħω − 1.86 meV: the quartic anharmonicity of −V_Γ(r), reproduced by")
println("first-order perturbation theory to 0.03 meV (see 03_selftest.jl Part E).")
println("The paper's 56.74 meV is the harmonic value (C24)/(C25).")
