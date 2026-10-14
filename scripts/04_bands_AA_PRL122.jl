# =============================================================================
# 04_bands_AA_PRL122.jl -- AA-stacked tMoTe2 with the parameters of
#                          PRL 122, 086402 (2019), first harmonics, at three angles.
#
#     julia scripts/04_bands_AA_PRL122.jl
#
# writes figs/prl122_AA_1.2deg.pdf    -- the Letter's own Fig. 3(a) angle
#        figs/prl122_AA_1.765deg.pdf  -- where its theta_1* gap closing happens
#        figs/prl122_AA_3.89deg.pdf   -- PRB 109's fitting angle, for comparison
#
# Same route as scripts/01_bands_AA.jl (Gamma_M -> K_M -> M_M -> Gamma_M), same
# energy zero (the +-K valence band maximum) and the same panel layout.  What
# changes is the angle and where the numbers come from: PRL 122 fits the
# AA-aligned bilayer, and draws its moire bands at theta = 1.2 deg.
#
# This script draws the figure; it is *not* the correctness test (that is
# scripts/03_selftest.jl, which runs at 3.89 deg).
# =============================================================================

using Printf
using Plots
gr()

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD

const NSHELL = 5
const FIGDIR = joinpath(ROOT, "figs")

# A choice of *this figure*, not of the model.
const C_K = RGB(0.80, 0.12, 0.12)     # ±K valley, as in 01/02

# PRL 122, 086402 (2019), just below Eq. (3): "(V, psi, w) ~ (8 meV, -89.6 deg,
# -8.5 meV) for MoTe2".  Two notes on provenance:
#
#   * the Letter defers m* to its Supplemental Material, so 0.62 (the monolayer
#     MoTe2 K-valley hole mass) is taken from outside it;
#   * V2 = w2 = 0 is the point of "FH": the Letter keeps only the first
#     harmonics, and PRB 109 confirms that this fit carries the top two valence
#     bands in each valley and no more.
#
# These numbers are a property of the *AA-aligned* bilayer -- the Letter fits
# them to the three high-symmetry interlayer displacements -- so they carry no
# preferred twist angle, and one payload is reused for both panels below.  Each
# panel still gets its own ParamRef naming the geometry it claims, so
# `check_parameters` sees a record that matches the geometry in hand.
const PRL122_PARAMS = KParams(mstar = 0.62, V = 8.0, psi = deg2rad(-89.6), w = -8.5,
                              V2 = 0.0, w2 = 0.0, label = "AA +/-K FH (PRL 122)")

const PRL122_SOURCE = "PRL 122, 086402 (2019): (V, psi, w) fitted to the AA-aligned " *
                      "bilayer, its moire bands drawn at theta = 1.2 deg.  m* is quoted " *
                      "only in its Supplemental Material"

# 1.2 deg is the Letter's Fig. 3(a).  1.765 deg is where *this* parameter set
# closes the gap between its 2nd and 3rd bands at Γ_M; the Letter quotes
# theta_1* ~ 1.74 deg for that closing, so those two panels are the same physics
# on either side of it.
#
# 3.89 deg is not an angle of PRL 122 at all -- it is PRB 109's, the angle at
# which the FH+SH fit of scripts/01_bands_AA.jl is defined.  Drawing this same
# PRL 122 payload there is what makes the two parameter sets comparable band by
# band, and it is also the angle scripts/03_selftest.jl runs at.
#
# Each panel carries its own window, because the energy scale is not shared:
# at 1.2 deg the top two ±K bands span a few meV and the window is a zoom, as in
# the Letter's Fig. 3(a); at 3.89 deg the three bands span 65 meV.
const PANELS = (("prl122_AA_1.2deg",    1.2,   (-14.0, 1.0)),
                ("prl122_AA_1.765deg",  1.765, (-14.0, 1.0)),
                ("prl122_AA_3.89deg",   3.89,  (-70.0, 4.0)))

for (tag, theta_deg, ylim) in PANELS
    prl122 = ParamSet(
        ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = theta_deg,
                 stacking = :AA, valley = :K, harmonics = :FH, provenance = FITTED,
                 source = PRL122_SOURCE),
        PRL122_PARAMS)

    # every argument is written out: this script doubles as the usage example,
    # and the defaults of `moire_geometry` are the same numbers (A0_MOTE2 ==
    # A0_DEFAULT)
    geo = moire_geometry(a0 = A0_MOTE2, theta_deg = theta_deg, Nshell = NSHELL)
    # the route of 01: Γ_M -> K_M -> M_M -> Γ_M; one k-set for the energy zero
    # and for every curve
    path = bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])
    basK = build_basis(geo, :AA_K)
    shift = vbm(geo, basK, prl122, :AA_K, path)   # energy zero = ±K valence band max
    # three bands, not two: band 3 is the one band 2 closes onto at theta_1*,
    # so it has to be on the plot for the lower panel to say anything.  Band 4
    # falls below every window here (it reaches -81 meV at 3.89 deg).
    cK = valence_bands(geo, basK, prl122, :AA_K, path; nbands = 3)

    # Band 1 is the topmost, as PRL 122 numbers them.
    @printf("%s: a_M = %.2f Å,  basis dimension %d\n", tag, geo.aM, length(basK.Qs))
    for j in size(cK.E, 2):-1:1
        @printf("  band %d: bandwidth %.4f meV,  at Γ_M %.4f,  at K_M %.4f,  at M_M %.4f\n",
                size(cK.E, 2) - j + 1, maximum(cK.E[:, j]) - minimum(cK.E[:, j]),
                cK.E[1, j] - shift, cK.E[cK.nodes[2], j] - shift, cK.E[cK.nodes[3], j] - shift)
    end
    for j in size(cK.E, 2):-1:2
        @printf("  gap between bands %d and %d: min over the path %.4f meV (at Γ_M %.4f)\n",
                size(cK.E, 2) - j + 1, size(cK.E, 2) - j + 2,
                minimum(cK.E[:, j] .- cK.E[:, j-1]), cK.E[1, j] - cK.E[1, j-1])
    end

    # --- the panel ---------------------------------------------------------
    # `$...$` switches on GR's LaTeX subset (GR's own parser, not real LaTeX):
    # `\Gamma_M` gives Γ with a proper subscript M, while a bare "Γ_M" would
    # print the underscore literally.  `raw"..."` is needed so the backslash
    # survives Julia's own string escaping -- and it does *not* interpolate, so
    # the angle is appended with `*` rather than written as `$theta_deg` inside
    # a raw string.  (A stray `$` also breaks GR's math-mode pairing: the two
    # dollars of `MoTe$_2$` must be the only ones in the raw part.)
    ticks = (cK.d[cK.nodes], [raw"$\Gamma_M$", raw"$K_M$", raw"$M_M$", raw"$\Gamma_M$"])
    plt = plot(; ylabel = raw"$E - E_{VBM}$  (meV)", ylim = ylim,
               title = raw"AA-stacked MoTe$_2$, FH (PRL 122), θ = " * "$(theta_deg)°",
               titlefontsize = 8, legend = :bottomright, legendfontsize = 8,
               titlelocation = (0.5, 1.02),
               size = (640, 420), xticks = ticks, framestyle = :box, grid = false,
               left_margin = 8Plots.mm)
    for j in 1:size(cK.E, 2)
        plot!(plt, cK.d, cK.E[:, j] .- shift; label = j == 1 ? "±K valley" : "",
              color = C_K, linestyle = :solid, linewidth = 1.4)
    end
    vline!(plt, cK.d[cK.nodes]; color = :gray70, linewidth = 0.6, label = "")
    hline!(plt, [0.0]; color = :gray70, linewidth = 0.6, label = "")
    mkpath(FIGDIR)
    savefig(plt, joinpath(FIGDIR, tag) * ".pdf")
    # savefig(plt, joinpath(FIGDIR, tag) * ".png")   # PDF only, by request
end

println("→ figs/prl122_AA_1.2deg.pdf, figs/prl122_AA_1.765deg.pdf, figs/prl122_AA_3.89deg.pdf")
println()
println("There is no Gamma-valley curve here: PRL 122 has no Gamma parameters, its")
println("model is the ±K two-band one.")
println()
println("At 1.2 deg the moire bands are flat on the scale of the 3.89 deg figures --")
println("the two topmost ±K bands span a few meV -- so those two panels share a zoom")
println("window, as in PRL 122's own Fig. 3(a), and their band 1-2 gaps can be read")
println("off directly.  The 3.89 deg panel needs a window six times taller, because")
println("that is where this parameter set actually spreads out.")
println()
println("Only the top two ±K bands are what the FH fit is meant to carry (PRB 109);")
println("band 3 is drawn because it is the one band 2 closes onto at Γ_M, the closing")
println("the Letter calls theta_1*.  To see that 2-3 closing itself, which is 0.007")
println("meV wide at 1.765 deg, the window has to be zoomed to about 1 meV.  At 3.89")
println("deg the same two bands are 25 and 13 meV wide and the closing is long gone.")
