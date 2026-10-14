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
using Plots
gr()

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD

const NSHELL = 5          # converged for the ±K bands (see 03_selftest.jl Part D)
const FIGDIR = joinpath(ROOT, "figs")

# Fig. 4 colours.  A choice of *this figure*, not of the model, so it is written
# here rather than in src/.
const C_K = RGB(0.80, 0.12, 0.12)     # ±K valley   (red in Fig. 4)
const C_G = RGB(0.95, 0.55, 0.05)     # Γ valley    (orange in Fig. 4)

# The paper fits 2 valence bands (one per valley, each doubly degenerate) with
# the FH parameters and 3 with FH+SH (App. B 3), so that is what is drawn.
for (tag, prmK, nbandK, title, leg) in
        (("fig4a_AA_FH",   AA_K_FH,   2, raw"AA-stacked MoTe$_2$, FH, θ = 3.89°",    :bottomright),
         ("fig4b_AA_FHSH", AA_K_FHSH, 3, raw"AA-stacked MoTe$_2$, FH+SH, θ = 3.89°", :bottomright))
    # every argument is written out: this script doubles as the usage example,
    # and `moire_geometry`'s defaults are the same numbers (A0_MOTE2 ==
    # A0_DEFAULT, THETA_DEFAULT_DEG) -- spelled here so the geometry is visible
    # rather than implied.
    geo = moire_geometry(a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG, Nshell = NSHELL)
    # Fig. 4's route Γ_M -> K_M -> M_M -> Γ_M; one k-set for the energy zero
    # and for every curve
    path = bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])
    basK = build_basis(geo, :AA_K)
    shift = vbm(geo, basK, prmK, :AA_K, path)      # energy zero = ±K valence band max
    cK = valence_bands(geo, basK, prmK, :AA_K, path; nbands = nbandK)
    cG = valence_bands(geo, build_basis(geo, :AA_Gamma), AA_GAMMA, :AA_Gamma, path;
                       nbands = 6)                 # the Γ basis carries spin
    @printf("%s: ±K VBM = %.3f meV,  top Γ band = %.3f meV (relative to the VBM)\n",
            tag, shift, maximum(cG.E) - shift)

    # --- the panel --------------------------------------------------------
    # `$...$` switches on GR's LaTeX subset (it is GR's own parser, not real
    # LaTeX): `\Gamma_M` gives Γ with a proper subscript M, while a bare
    # "Γ_M" would print the underscore literally.  `raw"..."` is needed so the
    # backslash survives Julia's own string escaping.
    ticks = (cK.d[cK.nodes], [raw"$\Gamma_M$", raw"$K_M$", raw"$M_M$", raw"$\Gamma_M$"])
    plt = plot(; ylabel = raw"$E - E_{VBM}$  (meV)", ylim = (-80.0, 5.0), title = title,
               titlefontsize = 8, legend = leg, legendfontsize = 8,
               titlelocation = (0.5, 1.02),
               size = (640, 420), xticks = ticks, framestyle = :box, grid = false,
               left_margin = 8Plots.mm)
    for (lab, col, ls, dat) in (("±K valley", C_K, :solid, cK),
                                ("Γ valley (unfitted)", C_G, :dash, cG))
        for j in 1:size(dat.E, 2)
            plot!(plt, dat.d, dat.E[:, j] .- shift; label = j == 1 ? lab : "",
                  color = col, linestyle = ls, linewidth = 1.4)
        end
    end
    vline!(plt, cK.d[cK.nodes]; color = :gray70, linewidth = 0.6, label = "")
    hline!(plt, [0.0]; color = :gray70, linewidth = 0.6, label = "")
    mkpath(FIGDIR)
    savefig(plt, joinpath(FIGDIR, tag) * ".pdf")
    # savefig(plt, joinpath(FIGDIR, tag) * ".png")   # PDF only, by request
end

println("→ figs/fig4a_AA_FH.pdf, figs/fig4b_AA_FHSH.pdf")
println()
println("Note: the Γ-valley parameters for AA stacking are not fitted anywhere in")
println("the paper (App. B 3 fits the ±K model alone, because for θ ≤ 5.09° the Γ")
println("valley lies below the 6th valence band).  AA_GAMMA reuses the AB values")
println("and its offset is a documented placeholder, so the orange curves sit well")
println("below the plotted window — which is exactly the paper's statement.")
