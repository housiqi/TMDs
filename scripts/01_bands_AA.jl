# =============================================================================
# 01_bands_AA.jl -- Fig. 4(a) and 4(b) of PRB 109, 205121 (2024): AA-stacked
#                   tMoTe2 at theta = 3.89 deg, first harmonics and
#                   first + second harmonics.
#
#     julia scripts/01_bands_AA.jl
#
# writes figs/fig4a_AA_FH.{pdf,png} and figs/fig4b_AA_FHSH.{pdf,png}
#
# The ±K bands are red, the Γ-valley bands orange, the energy zero is the ±K
# valence band maximum, and the path is Γ_M -> K_M -> M_M -> Γ_M -- the layout
# of the paper's Fig. 4.  This script draws the figure; it is *not* the
# correctness test (that is scripts/03_selftest.jl).
# =============================================================================

using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD
include(joinpath(@__DIR__, "plot_common.jl"))

const NSHELL = 5          # converged for the ±K bands (see 03_selftest.jl Part D)
const FIGDIR = joinpath(ROOT, "figs")

# The paper fits 2 valence bands (one per valley, each doubly degenerate) with
# the FH parameters and 3 with FH+SH (App. B 3), so that is what is drawn.
for (tag, prmK, nbandK, title) in
        (("fig4a_AA_FH",   AA_K_FH,   2, "AA stacking, FH,   θ = 3.89°"),
         ("fig4b_AA_FHSH", AA_K_FHSH, 3, "AA stacking, FH+SH, θ = 3.89°"))
    geo = moire_geometry(theta_deg = THETA_DEFAULT_DEG, Nshell = NSHELL)
    shift = vbm(geo, :AA_K, prmK)        # energy zero = ±K valence band maximum
    cK = valence_bands(geo, :AA_K, prmK; nbands = nbandK)
    cG = valence_bands(geo, :AA_Gamma, AA_GAMMA; nbands = 6)   # spinful basis
    @printf("%s: ±K VBM = %.3f meV,  top Γ band = %.3f meV (relative to the VBM)\n",
            tag, shift, maximum(cG.E) - shift)
    fig4_panel([("±K valley", C_K, :solid, cK, shift),
                ("Γ valley (unfitted)", C_G, :dash, cG, shift)],
               joinpath(FIGDIR, tag); title = title, ylim = (-80.0, 5.0))
end

println("→ figs/fig4a_AA_FH.pdf(.png), figs/fig4b_AA_FHSH.pdf(.png)")
println()
println("Note: the Γ-valley parameters for AA stacking are not fitted anywhere in")
println("the paper (App. B 3 fits the ±K model alone, because for θ ≤ 5.09° the Γ")
println("valley lies below the 6th valence band).  AA_GAMMA reuses the AB values")
println("and its offset is a documented placeholder, so the orange curves sit well")
println("below the plotted window — which is exactly the paper's statement.")
