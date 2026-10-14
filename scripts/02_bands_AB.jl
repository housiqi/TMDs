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
using Plots
gr()

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD

const NSHELL = 7          # the Γ-valley oscillator needs Nshell ≥ 7; it is far
const FIGDIR = joinpath(ROOT, "figs")   # more localized than the ±K states

# Fig. 4 colours -- a choice of this figure, not of the model.
const C_K = RGB(0.80, 0.12, 0.12)     # ±K valley   (red in Fig. 4)
const C_G = RGB(0.95, 0.55, 0.05)     # Γ valley    (orange in Fig. 4)

# every argument is written out: this script doubles as the usage example, and
# the defaults of `moire_geometry` happen to be the same numbers
geo = moire_geometry(a0 = A0_MOTE2, theta_deg = THETA_DEFAULT_DEG, Nshell = NSHELL)
# Fig. 4's route Γ_M -> K_M -> M_M -> Γ_M; one k-set for the energy zero and
# every curve
path = bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])
basK = build_basis(geo, :AB_K)
basG = build_basis(geo, :AB_Gamma)
shift = vbm(geo, basK, AB_K, :AB_K, path)
cK = valence_bands(geo, basK, AB_K, :AB_K, path; nbands = 4)
# the Γ-valley basis carries spin explicitly, so the top 6 states are the
# (0,0) doublet plus the (1,0),(0,1) quartet -- the paper's 2-fold + 4-fold
cG = valence_bands(geo, basG, AB_GAMMA, :AB_Gamma, path; nbands = 6)
@printf("±K VBM = %.3f meV;  Γ flat bands at %.3f and %.3f meV below it\n",
        shift, maximum(cG.E) - shift, minimum(cG.E) - shift)
@printf("ħω (C25) = %.3f meV;   exact lattice splitting = %.3f meV\n",
        hbar_omega(geo, AB_GAMMA), maximum(cG.E) - minimum(cG.E))

# --- panel 1: the whole valence window, both valleys ------------------------
# `$...$` switches on GR's LaTeX subset (GR's own parser, not real LaTeX):
# `\Gamma_M` gives Γ with a proper subscript M, while a bare "Γ_M" would print
# the underscore literally.  `raw"..."` keeps the backslash from Julia.
ticks = (cK.d[cK.nodes], [raw"$\Gamma_M$", raw"$K_M$", raw"$M_M$", raw"$\Gamma_M$"])
plt = plot(; ylabel = raw"$E - E_{VBM}$  (meV)",
           ylim = (-100.0, 5.0), title = raw"AB-stacked MoTe$_2$, θ = 3.89°",
           titlefontsize = 8, legend = :bottomright, legendfontsize = 8,
           titlelocation = (0.5, 1.02),
           size = (640, 420), xticks = ticks,
           framestyle = :box, grid = false, left_margin = 8Plots.mm)
for (lab, col, dat) in (("±K valley", C_K, cK), ("Γ valley", C_G, cG))
    for j in 1:size(dat.E, 2)
        plot!(plt, dat.d, dat.E[:, j] .- shift; label = j == 1 ? lab : "",
              color = col, linestyle = :solid, linewidth = 1.4)
    end
end
vline!(plt, cK.d[cK.nodes]; color = :gray70, linewidth = 0.6, label = "")
hline!(plt, [0.0]; color = :gray70, linewidth = 0.6, label = "")
mkpath(FIGDIR)
savefig(plt, joinpath(FIGDIR, "fig4c_AB") * ".pdf")
# savefig(plt, joinpath(FIGDIR, "fig4c_AB") * ".png")   # PDF only, by request

# --- panel 2: the two flat bands alone, on the paper's −30 / −86 meV scale ---
# they carry the whole Γ-valley story
plt2 = plot(; ylabel = raw"$E - E_{VBM}$  (meV)",
            ylim = (-95.0, -20.0),
            # raw"..." cannot interpolate $NSHELL, hence the concat: the LaTeX part
            # stays raw, the counter is interpolated.
            title = raw"AB-stacked MoTe$_2$" * ": Γ-valley flat bands (Nshell = $NSHELL)",
            titlefontsize = 8, legend = :bottomright, legendfontsize = 8,
            titlelocation = (0.5, 1.02),
            size = (640, 420), xticks = ticks, framestyle = :box, grid = false,
            left_margin = 8Plots.mm)
for j in 1:size(cG.E, 2)
    plot!(plt2, cG.d, cG.E[:, j] .- shift; label = j == 1 ? "Γ valley flat bands" : "",
          color = C_G, linestyle = :solid, linewidth = 1.4)
end
vline!(plt2, cG.d[cG.nodes]; color = :gray70, linewidth = 0.6, label = "")
hline!(plt2, [0.0]; color = :gray70, linewidth = 0.6, label = "")
savefig(plt2, joinpath(FIGDIR, "fig4c_AB_gamma_zoom") * ".pdf")
# savefig(plt2, joinpath(FIGDIR, "fig4c_AB_gamma_zoom") * ".png")   # PDF only, by request

println("→ figs/fig4c_AB.pdf, figs/fig4c_AB_gamma_zoom.pdf")
println()
println("The two flat bands are the (n_x,n_y) = (0,0) and (1,0),(0,1) levels of the")
println("harmonic oscillator of App. C 3 — single and double spinless C3 = e^{±i2π/3},")
println("i.e. 2-fold and 4-fold with spin.  Their exact splitting is 54.84 meV, which")
println("is ħω − 1.86 meV: the quartic anharmonicity of −V_Γ(r), reproduced by")
println("first-order perturbation theory to 0.03 meV (see 03_selftest.jl Part E).")
println("The paper's 56.74 meV is the harmonic value (C24)/(C25).")
