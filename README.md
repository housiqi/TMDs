# MoireTMD — continuum models of twisted bilayer systems

A Julia package for the single-particle continuum models of
> Phys. Rev. B **109**, 205121 (2024).

All four models of that paper are implemented: **AA stacking** (App. B) and
**AB stacking** (App. C), each with the **±K valley** and the **Γ valley**. Every model is diagonalised in the plane-wave basis.

---

## band calculation

Every script for band calculation has the same steps.  This is the entire interface:

```julia
using Printf, Plots; gr()

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MoireTMD.jl"))
using .MoireTMD

# 1. geometry -- a0 and theta_deg are pure inputs; nothing here is fitted
geo = moire_geometry(a0 = A0_MOTE2, theta_deg = 3.89, Nshell = 5)

# 2. the k-route; one k-set, shared by the energy zero and every curve
path = bz_path(geo, [[0.0, 0.0], geo.K_M, geo.Mp_M, [0.0, 0.0]]; nseg = [40, 20, 40])

# 3. the basis belongs to the model, not to the parameters
basK = build_basis(geo, :AA_K)

# 4. the parameters: a payload, plus a record of where the numbers came from
mine = ParamSet(
    ParamRef(material = :MoTe2, a0 = A0_MOTE2, theta_deg = 3.89, stacking = :AA,
             valley = :K, harmonics = :FH, provenance = EXPLORATORY,
             source = "my own numbers"),
    KParams(mstar = 0.62, V = 8.0, psi = deg2rad(-89.6), w = -8.5,
            V2 = 0.0, w2 = 0.0, label = "AA ±K FH (mine)"))

# 5. energy zero = the ±K valence band maximum, on this same k-set
shift = vbm(geo, basK, mine, :AA_K, path)

# 6. the top `nbands` valence bands; nothing is shifted yet
cK = valence_bands(geo, basK, mine, :AA_K, path; nbands = 3)

# 7. plot: `d` is the cumulative path length, `E` is nk × nbands
plt = plot(; ylabel = raw"$E - E_{VBM}$  (meV)", ylim = (-14.0, 1.0),
           legend = :bottomright)
for j in 1:size(cK.E, 2)
    plot!(plt, cK.d, cK.E[:, j] .- shift; label = j == 1 ? "±K valley" : "")
end
vline!(plt, cK.d[cK.nodes]; color = :gray70, linewidth = 0.6, label = "")
savefig(plt, "mine.pdf")
```

Things worth knowing before you change it:

* **`model`** is `:AA_K`, `:AB_K`, `:AA_Gamma` or `:AB_Gamma`, and
  `build_basis(geo, model)` picks the matching basis.  The Γ basis carries spin
  explicitly (4 states per `G_M`), so `nbands = 6` there is the top doublet plus
  the `(1,0),(0,1)` quartet — see 02.
* **`Nshell`** is the truncation.  ±K converges by 5–6; the Γ valley needs
  ≥ 7 (see *Notes*).  Every `valence_bands` / `vbm` call is a *sampling*, so
  both take the k-set from you and have no default.
* **`nbands`** counts from the top (eigenvalues ascend), so it selects the
  topmost valence bands — which is what these fits are for.  Ask for more than
  the basis has and you get an `ArgumentError`, not a silent pad.
* **`raw"..."`** switches on GR's LaTeX subset (`$...$` is math mode) and keeps
  the backslash from Julia.  `raw"$\Gamma_M$"` gives Γ with a subscript;
  a bare `"Γ_M"` prints the underscore.  `raw` does *not* interpolate —
  concatenate with `*` to mix in a computed value.

---

## `ParamSet` — the parameters, and where they came from

Table I and Table II are **fits to DFT at one material and one twist angle**
(MoTe₂ at θ = 3.89°), so a parameter set is never just a bag of numbers.  Each
one is a `ParamSet`: a physics payload wrapped with a `ParamRef` recording its
origin and the geometry it is valid at.

```julia
KParams(; mstar, V, psi, w, V2 = 0.0, w2 = 0.0, label = "")       # ±K valleys
GammaParams(; mstar, V, psi, w, E0 = 0.0, label = "")             # Γ valley

ParamRef(; material, a0, theta_deg, stacking, valley, harmonics,
         provenance, source, calibrated_at)                       # where it came from
provenance ∈ (FITTED, CALIBRATED, PLACEHOLDER, EXPLORATORY)
```

`payload(ps)` and `paramref(ps)` unwrap the two halves; every public function
takes either a `ParamSet` or a bare payload, so a quick hand-built `KParams`
still works.

| set | parameters | provenance | valid at |
|---|---|---|---|
| `AA_K_FH`   | m\* = 0.60, V = 16.5, ψ = −105.9°, w = −18.8 | FITTED, Table I | MoTe₂, 3.89° |
| `AA_K_FHSH` | m\* = 0.62, V = 7.94, ψ = −88.43°, w = −10.77, V₂ = 20.00, w₂ = 10.21 | FITTED, Table I | MoTe₂, 3.89° |
| `AB_K`      | m\* = 0.62, V = 53, ψ = −56°, w = 0 | FITTED, Table II | MoTe₂, 3.89° |
| `AB_GAMMA`  | m_Γ\* = 10, V_Γ = 72, ψ_Γ = 0°, w_Γ = 300, E_Γ = −640.13 | CALIBRATED — Table II except E_Γ | MoTe₂, 3.89° |
| `AA_GAMMA`  | *(no fitted values exist — see below)* | PLACEHOLDER | MoTe₂, 3.89° |

```julia
p   = paper_fit(:AA, :K, :FHSH)          # keys are (stacking, valley, harmonics)
geo = moire_geometry(MOTE2; theta_deg = 3.89)
check_parameters(geo, p; model = :AA_K)  # throws if p is not valid at this geometry
```

**The guard.**  `check_parameters` runs by itself at the top of `hamiltonian`
and `transcribe_hamiltonian` — the two places where a parameter set meets a
geometry — and raises if `a0`, `theta_deg`, or the stacking/valley it was fitted
for disagrees with what is being computed.  The reason is in the paper's own
Table V: the DFT top-valence bandwidth moves from 12.8 meV at 3.89° to 133.3 meV
at 9.43°, and the FH+SH second harmonics V₂, w₂ are exactly the terms standing
in for the relaxation behind that growth.  Carried to another angle, a fitted set
produces smooth, plausible, *wrong* bands — and since `moire_geometry` takes `a0`
and `theta_deg` as free inputs, that used to happen silently.

Two escape hatches, for when crossing the boundary is deliberate:

* `exploratory(ps)` — the same payload with the geometry claim withdrawn, per
  record (that is what the example above does);
* `enforce_param_checks!(false)` — the same thing session-wide.

A record that makes no claim (`a0` or `theta_deg` = `NaN`, or provenance
`EXPLORATORY`) is skipped, so your own `KParams` is unconstrained by default.

**Adding a material** — WSe₂, say — is a `Material` for its `a0` plus one
`ParamSet` per model carrying *its own* V, ψ, w, … from whatever source you are
using; there is no registry to edit.  PRB 109 fits MoTe₂ only, so nothing here
carries over except the geometry.

**Two places where the paper is silent**, both flagged in `parameters.jl`:

1. **AA Γ valley has no fitted parameters at all.**  App. B 3 fits the ±K-valley
   model alone, because for θ ≤ 5.09° the Γ-valley valence bands lie below the
   6th valence band.  `AA_GAMMA` reuses the AB values (at ψ_Γ = 0 the (B38) and
   (C13) potentials are the same expression) and its offset is a **placeholder**
   160 meV below the ±K VBM.  They are not fits.
2. **`E_Γ` is tabulated nowhere.**  It only slides the whole Γ-valley manifold
   relative to the ±K valleys (App. B 3), so it is a free parameter.  The stored
   value is *calibrated*, not fitted, from the paper's statement that the first
   ultra-flat Γ band lies 30 meV below the ±K VBM; `E0_relative_to_vbm`
   recomputes it, and the self-test reproduces the stored number to 3×10⁻⁴ meV.
   Only the *position* depends on it, never the band structure or the flat-band
   splitting.

---

## Files

```
Project.toml                 deps: LinearAlgebra, Printf, Plots
src/MoireTMD.jl              module entry point
src/geometry.jl              a_0, θ → a_M, b_{M,i}, g_i, g_{2i}, q_i, q_{2i}, Q sets, k-path
src/parameters.jl            ParamSet / ParamRef records, Table I / Table II constants,
                             the paper_fit registry, and the geometry-validity guard
src/potentials.jl            analytic Fourier components of (B21)/(B38)/(C4)/(C13)
src/kmodels.jl               assembles h_QQ'(k): AA-K / AB-K / AA-Γ / AB-Γ
src/transcribe.jl            verbatim transcription of the printed (B27)/(C8)/(C18)
src/symmetry.jl              every symmetry operator as an explicit basis matrix, + the suite
src/bands.jl                 diagonalize_point, vbm, valence_bands, the E_Γ calibration rule
scripts/01_bands_AA.jl       Fig. 4(a)/(b)
scripts/02_bands_AB.jl       Fig. 4(c) + Γ-valley zoom
scripts/03_selftest.jl       the whole self-test
scripts/04_bands_AA_PRL122.jl  PRL 122, 086402 (2019) parameters
```
---

## The four models

```
h_QQ'(k)  =  δ_QQ' [ −ħ²(k−Q)²/(2m*) + (ε/2)·σ(Q) ]        kinetic + displacement field
          +  intralayer potential  V(r)                     (B21) / (B38) / (C4) / (C13)
          +  interlayer tunnelling t(r)                     (B21) / (B38) / (B4) / (C11)
```

| model | stacking | valley | real-space source | printed momentum form |
|---|---|---|---|---|
| `:AA_K`     | AA | ±K    | (B15), (B17), (B20), (B21) | (B27), Eq. (3) |
| `:AB_K`     | AB | ±K    | (C4), Eq. (4)/(5)          | (C8), Eq. (5)  |
| `:AA_Gamma` | AA | Γ     | (B28), (B32), (B38)        | — *(not printed)* |
| `:AB_Gamma` | AB | Γ     | (C11), (C13)               | (C18), Eq. (7) |

**Matrix-element rule.**  (B22)/(C16) expand the *field operator*, not the
state, so a state created by `c†_{η,l,k−Q}` carries the conjugate kernel
`e^{+i(k−Q)·r}`, and for `O(r) = Σ_G O_G e^{−ηiG·r}`

```
⟨k,Q|O|k,Q'⟩ = ∫dr e^{−i(k−Q)·r} O(r) e^{+i(k−Q')·r} = O_{η, η(Q−Q')} ,
```

which is exactly the index ordering the printed (B27)/(C8)/(C18) use.  Reading
the (B22) kernel as the wavefunction instead gives `O_{η,η(Q'−Q)}` — the plain
transpose, i.e. the same model with ψ → −ψ.  Both are implemented
(`convention = :paper` / `:transposed`) and, `H` being Hermitian, they have
identical spectra.  `:paper` is the default.

---

## How correctness is judged

The criteria are the paper's own identities
and the appendices' own analytic results, all in `scripts/03_selftest.jl`.

* **A — operator identities** (residual ‖LHS−RHS‖ at a generic k, all passing):
  Hermiticity; Bloch periodicity `H(k + b_{M,i}) = H(k)`; C₃ covariance and the
  Γ_M representations {1, ω, ω²}; effective inversion (B18)/(B39); the mirror-T
  composites (B10)/(B14)/(C15); valley TR and the AB effective TR (C5)/(C14);
  `H_{+1}(k) = H_{−1}(k)` identically for AB±K; exact 2-fold spin degeneracy at
  Γ.  Two are *predicted failures* and are reported as such rather than
  asserted: effective inversion for AB±K, and the AB effective 𝒦 for AA.
* **B — main source ↔ cross source.**  `kmodels.jl` builds `h_QQ'(k)` from the
  analytic Fourier components of the real-space potentials; `transcribe.jl`
  re-derives the same matrix by transcribing the printed (B27)/(C8)/(C18)
  verbatim, scanning every basis pair, never touching the Fourier-table
  machinery.  They agree to `0.00e+00` for all five parameter sets at both
  valleys.
* **C — quantitative.**  η = ±1 have identical spectra point by point for AA±K;
  the `w = 0` layer degeneracy is fourfold along Γ_M–M_M but only twofold at K_M
  (the paper's Fig. 4(c) pattern); the top 6 ±K bands move by ≤ 3×10⁻⁵ meV
  between Nshell = 4, 5, 6; `ħω` from (C25) reproduces 56.76 meV to 0.0014 meV.
* **D — the Γ-valley flat bands** against App. C 3 (next section).
* **E — the parameter records** (Part G, which runs after A–F).  Every fitted
  number and the whole band structure are frozen — three Σ|E| checksums along
  Γ_M–K_M–M_M–Γ_M, the ±K VBM, `ħω`, the anharmonic shift — and asserted
  *bit-identical* to the values from before the `ParamSet` refactor, because
  provenance metadata has to be physically inert.  It then exercises the
  registry, each of the three mismatch failures, both escape hatches, and the
  re-stamping of the `E_Γ` calibration.
* **ε ≠ 0 is the lever.**  Every identity whose operator exchanges the layers
  **must** break at ε ≠ 0, since (B19)/(B40)/(C15) make the displacement field
  `(−1)^l ε/2` flip sign under the exchange.  Part C drives them all at
  ε = 10 meV and confirms they do.  At ε = 0 the Γ-valley models are
  layer-exchange symmetric *on their own*, so a layer-flip error in a symmetry
  permutation would pass either way.

---

## The Γ-valley flat bands, and the paper's 56.74 meV

App. C 3 reduces the AB Γ-valley model to decoupled harmonic oscillators:
`m_Γ* = 10 m_e` gives an oscillator length ≈ 3.66 Å, tiny inside the 52 Å moiré
cell, so the top states are the `(n_x, n_y)` levels of `−V_Γ(r)` about the AA
site, with `ħω = ħ√(3V_Γ|b_{M,1}|²/m_Γ*) = 56.76 meV` (C25).  Then
`(0,0)` is one spinless state with C₃ = 1 (2-fold with spin) and
`(1,0),(0,1)` are two with C₃ = e^{±i2π/3} (4-fold) — exactly what the
diagonalisation returns, at −30.00 and −84.90 meV, with bandwidths 0.003 and
0.011 meV over the whole MBZ.

**The splitting is 54.84 meV, not 56.74.**  The difference is the quartic
anharmonicity (C23) drops: expanding (C20) exactly,
`−V_Γ(r) = −6V_Γ + (3/2)V_Γ|b|²r² − (3/32)V_Γ|b|⁴r⁴ + …`, and first-order
perturbation theory with the 2D oscillator moments gives

```
ΔE₁ − ΔE₀ = −(3/8) V_Γ |b|⁴ l⁴ = −(1/4) ħ²|b_{M,1}|²/(2m_Γ*) = −1.865 meV,
```

independent of `V_Γ`.  So `56.761 − 1.865 = 54.897` meV, against 54.863 from the
numerical diagonalisation at Nshell = 7 — reproduced to 0.034 meV, which is
`anharmonic_shift` in `parameters.jl`.  The paper's 56.74 is (C25)'s harmonic
value, consistent with App. C 3 describing the whole analysis as an
approximation.

---

## Notes and caveats

* **The Γ-valley label carries no layer information.**  In ±K the label
  `Q = G_M + η(−1)^l q₁` *encodes* the layer, so `Q → −Q` flips it automatically.
  The Γ label set is just `{G_M}`, so every operation that exchanges layers must
  do it in the permutation — (B39), (C14), (B14).
  `basis_perm(::GammaBasis, ...; layerswap = true)` is that switch.
* **Mirror × T acts on the momentum label as the *other* reflection.**  `C₂ᵧ𝒯`
  maps the label by C₂ₓ and `k → C₂ₓk`, while AB's `C₂ₓ𝒯` maps it by C₂ᵧ.  This is
  forced, not chosen: only that reflection takes `q₁ → −q₁`, which lies in the
  *other* family of (B23) and is therefore exactly the layer flip.
* **App. B and App. C disagree about which layer is `l = 0`** — a naming choice,
  equivalent to ε → −ε.  At ε = 0, where every band structure here is computed,
  it is irrelevant.  `layer_sign_kvalley` / `layer_sign_gamma` keep each
  appendix's own printed convention.
* **Truncation.**  Each set is cut in a way its own point group leaves invariant,
  so the cut never breaks a symmetry.  The Γ triangular set `{G_M}` uses the hex
  shell `|m1|, |m2|, |m1+m2| ≤ Nshell`.  The ±K set `{G_M ± q₁}` is *centred on
  ±q₁*, for which that same shell is **not** C₃-invariant, so it is cut by the
  C₃-and-inversion-invariant disk `|Q| ≤ Rcut` instead, with `Rcut` inherited
  from the hex shell so that `Nshell` still sets the size.  The Γ models need
  **Nshell ≥ 7**; their splitting converges slowly from above
  (62.23 / 56.53 / 55.09 / 54.86 / 54.84 / 54.84 meV for N = 4…9).
* **The *low* Γ-valley bands are not converged and are not shown.**  With the
  valence-band sign the Γ spectrum is unbounded below, and states below roughly
  −300 meV keep moving as Nshell grows.  Only the top of the Γ manifold is a
  converged statement; the plotted window lies entirely inside it.
* **No OCR error in (B27) or (C18).**  An earlier plan asserted two — that (B27)
  has the `±g_i` deltas of (B21) exchanged, and that (C18)'s sign differs from
  (C13).  Both are wrong: each is consistent with the real-space equations
  through the matrix-element rule above, which is also why `kmodels.jl` and
  `transcribe.jl` agree exactly.  The other OCR observations in that plan stand:
  Table II's columns are misaligned in the extracted text, and several equation
  numbers (B12, C22) are swallowed by the extraction.

---

## Reference

Jia, Zhang, Yu, Liu, Song, Wang, Yao, *Moiré fractional Chern insulators. I.
Origin of the insulating states and quantum anomalous Hall effect*,
Phys. Rev. B **109**, 205121 (2024).  Equation numbers in parentheses, e.g.
(B27), refer to its appendices B and C.
