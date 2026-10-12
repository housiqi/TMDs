# =============================================================================
# plot_common.jl -- shared plotting for 01_bands_AA.jl and 02_bands_AB.jl
#
# Reproduces the layout of the paper's Fig. 4: the high-symmetry path
# Γ_M -> K_M -> M_M -> Γ_M, red lines from the ±K valley and orange lines from
# the Γ valley, energies measured from the ±K valence band maximum.
# =============================================================================

using Plots
gr()

const C_K = RGB(0.80, 0.12, 0.12)     # ±K valley   (red in Fig. 4)
const C_G = RGB(0.95, 0.55, 0.05)     # Γ valley    (orange in Fig. 4)

"""
    valence_bands(geo, model, p; nbands, eta = 1, eps = 0.0)

Top `nbands` valence bands of `model` along the path, as a named tuple
`(d, nodes, E)` with `d` the cumulative path length and `E` an `nk × nbands`
matrix.  Nothing is shifted yet: the caller chooses the energy zero.
"""
function valence_bands(geo::MoireGeometry, model::Symbol, p; nbands::Int,
                       eta::Int = 1, eps::Real = 0.0)
    bd = diagonalize_path(geo, model, p; eta = eta, eps = eps)
    return (d = bd.kdist, nodes = bd.nodes, E = top_bands(bd, nbands))
end

"""
    fig4_panel(curves, fname; title, ylim, xlab)

Draw one Fig. 4-style panel.  `curves` is a vector of
`(label, color, linestyle, (d, nodes, E), shift)` tuples; every curve is plotted
as `E - shift` so that all of them share one energy zero.  Writes
`fname.pdf` and `fname.png`.
"""
function fig4_panel(curves, fname; title = "", ylim = (-80.0, 5.0),
                    xlab = "Γ_M        K_M        M_M        Γ_M")
    d, nodes = curves[1][4].d, curves[1][4].nodes
    ticks = (d[nodes], ["Γ_M", "K_M", "M_M", "Γ_M"])
    plt = plot(; xlabel = xlab, ylabel = "E − E_VBM  (meV)", ylim = ylim,
               title = title, titlefontsize = 10, legend = :topright,
               legendfontsize = 8, size = (640, 420), xticks = ticks,
               framestyle = :box, grid = false, left_margin = 8Plots.mm)
    for (lab, col, ls, dat, sh) in curves
        for j in 1:size(dat.E, 2)
            plot!(plt, dat.d, dat.E[:, j] .- sh; label = j == 1 ? lab : "",
                  color = col, linestyle = ls, linewidth = 1.4)
        end
    end
    vline!(plt, d[nodes]; color = :gray70, linewidth = 0.6, label = "")
    hline!(plt, [0.0]; color = :gray70, linewidth = 0.6, label = "")
    mkpath(dirname(fname))
    savefig(plt, fname * ".pdf")
    savefig(plt, fname * ".png")
    return plt
end
