# MoireTMD — continuum models of twisted bilayer MoTe₂

A Julia package reproducing the single-particle continuum models of

> **Jia, Zhang, Yu, Liu, Song, Wang, Yao, *Moiré fractional Chern insulators. I.
> Origin of the insulating states and quantum anomalous Hall effect*,
> Phys. Rev. B **109**, 205121 (2024)**

for twisted bilayer MoTe₂ at the twist angle θ = 3.89° used throughout that
paper.  All four models of the paper are implemented: **AA stacking** (App. B)
and **AB stacking** (App. C), each with the **±K valley** and the **Γ valley**.

The output is band structures in the style of the paper's Fig. 4.  Chern
numbers, quantum geometry, C₃ eigenvalues of individual bands and the FCI
analysis of the paper are explicitly *out of scope* here.

Everything is diagonalised in the plane-wave basis of (B22)/(C16),
`|k,Q⟩ ∝ exp(−i(k−Q)·r)` — exactly the route the paper takes.  There is no
real-space diagonalisation anywhere.

---

## Quick start

From the repository root, in the **default** Julia environment (so that the
already-installed `Plots` is visible):

```bash
julia scripts/03_selftest.jl     # correctness: symmetry identities + cross-checks
julia scripts/01_bands_AA.jl     # Fig. 4(a) and 4(b): AA stacking
julia scripts/02_bands_AB.jl     # Fig. 4(c): AB stacking
```

`03_selftest.jl` exits non-zero if any identity the paper predicts fails to
hold.  The figure scripts write `figs/*.pdf` and `figs/*.png`, which are the
deliverable to compare by eye — they are *not* the correctness test.

```julia
using MoireTMD
geo = moire_geometry(theta_deg = 3.89, Nshell = 7)
bas = build_basis(geo, :AB_Gamma)
E   = diagonalize_point(geo, bas, AB_GAMMA, :AB_Gamma, geo.K_M)
```

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
src/bands.jl                 diagonalize_point, vbm, and the E_Γ calibration rule
scripts/01_bands_AA.jl       Fig. 4(a)/(b)
scripts/02_bands_AB.jl       Fig. 4(c) + Γ-valley zoom
scripts/03_selftest.jl       the whole self-test
scripts/01,02 draw their own panels: the Fig. 4 route and colours are written
                             out in each script, not in src/ -- they belong to the figure
```

The two pre-existing example scripts in this repository
(`moireband_PRL121.jl`, `moireband_PRL122.jl`) are **untouched**.  (For the
record: both `include("parameters.jl")`, a file that does not exist in the
repository, and `moireband_PRL122.jl` uses `kappa1`/`kappa2`, which are never
defined anywhere — they do not run as they stand.)  Nothing here imports them.

---

## The four models, and where each matrix element comes from

For all four models the Hamiltonian is assembled as

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
state: `c†_{η,l,r} = V^{-1/2} Σ_k Σ_Q e^{−i(k−Q)·r} c†_{η,l,k−Q}`.  Both sides
are operators, so the state created by `c†_{η,l,k−Q}` carries the *conjugate*
kernel, `e^{+i(k−Q)·r}`, and for `O(r) = Σ_G O_G e^{−ηiG·r}`

```
⟨k,Q|O|k,Q'⟩ = ∫dr e^{−i(k−Q)·r} O(r) e^{+i(k−Q')·r} = O_{η, η(Q−Q')} ,
```

which is exactly the index ordering the printed (B26)/(B27) and (C8)/(C18) use.
Reading the kernel of (B22) as the wavefunction instead gives the opposite
ordering `O_{η,η(Q'−Q)}` — the plain transpose, i.e. the same model with ψ → −ψ.
Both readings are implemented (`convention = :paper` / `:transposed`); since `H`
is Hermitian, transpose = complex conjugate, so they have *identical spectra*,
and `‖H_paper − transpose(H_transposed)‖ = 0` exactly for every parameter set
(`03_selftest.jl`, Part B).  `:paper` is the default.

---

## Parameters: the fits, where they are valid, and the two places where the paper is silent

Table I and Table II are **fits to the DFT bands at one material and one twist
angle** — θ = 3.89°, MoTe₂ — and that is not a footnote.  The paper's own Table V
puts the DFT top-valence bandwidth at 12.8 meV at 3.89° and 133.3 meV at 9.43°,
and the FH+SH second harmonics V₂, w₂ are exactly the terms that stand in for the
lattice relaxation responsible for that growth.  So a fitted set carried to
another angle yields smooth, plausible, wrong bands.  The geometry side of this
package has no such restriction — `moire_geometry` takes `a0` and `theta_deg` as
pure inputs — which means the mismatch used to pass silently.  It no longer does.

Each set is a `ParamSet`: a physics payload (`KParams` / `GammaParams`) wrapped
with a `ParamRef` recording where the numbers came from and which geometry they
are valid at.

| set | parameters | provenance | valid at |
|---|---|---|---|
| `AA_K_FH`   | m\* = 0.60, V = 16.5, ψ = −105.9°, w = −18.8 | FITTED, Table I | MoTe₂, 3.89° |
| `AA_K_FHSH` | m\* = 0.62, V = 7.94, ψ = −88.43°, w = −10.77, V₂ = 20.00, w₂ = 10.21 | FITTED, Table I | MoTe₂, 3.89° |
| `AB_K`      | m\* = 0.62, V = 53, ψ = −56°, w = 0 | FITTED, Table II | MoTe₂, 3.89° |
| `AB_GAMMA`  | m_Γ\* = 10, V_Γ = 72, ψ_Γ = 0°, w_Γ = 300, E_Γ = −640.13 | CALIBRATED — Table II except E_Γ, see below | MoTe₂, 3.89° |
| `AA_GAMMA`  | *(no fitted values exist — see below)* | PLACEHOLDER | MoTe₂, 3.89° |

```julia
p   = paper_fit(:AA, :K, :FHSH)          # == AA_K_FHSH; keys are (stacking, valley, harmonics)
geo = moire_geometry(MOTE2; theta_deg = 3.89)
check_parameters(geo, p; model = :AA_K)  # throws if p is not valid at this geometry
E   = diagonalize_point(geo, build_basis(geo, :AA_K), p, :AA_K, [0.0, 0.0])
```

`check_parameters` also runs by itself at the top of `hamiltonian` and
`transcribe_hamiltonian` — the two points where a parameter set meets a geometry
— and raises if `a0`, `theta_deg`, or the stacking/valley the set was fitted for
disagrees with what is being computed (`03_selftest.jl` Part G asserts each of
those three failures).  Two escape hatches exist for when crossing the boundary
is deliberate:

* `exploratory(ps)` — the same payload with the geometry claim withdrawn, per
  record;
* `enforce_param_checks!(false)` — the same thing session-wide.

A record that makes no claim (`a0` or `theta_deg` = `NaN`, or provenance
`EXPLORATORY`) is skipped, so a hand-built `KParams` / `GammaParams` is
unconstrained by default.

**Adding a material** — WSe₂, say — is a `Material` for `a0` plus one `ParamSet`
per model carrying *its own* V, ψ, w, … from whatever source you are using; there
is no registry to edit.  PRB 109, 205121 (2024) fits MoTe₂ only, so nothing here
carries over except the geometry.

**Why these stay `const`.**  For an immutable record that is idiomatic and does
not restrict generality: any number of further `ParamSet`s can be defined
alongside them, and `const` only forbids rebinding the same *name*.  Dropping it
would be worse — a non-`const` global binding is typed `Any` at every use site.
What actually stood in the way of other materials and angles was not `const` but
the literals `3.52` and `3.89` sitting inside signature defaults, and those are
now `A0_DEFAULT` / `THETA_DEFAULT_DEG`.

**1. AA Γ valley has no fitted parameters at all.**  App. B 3 states that the
fits used the ±K-valley model alone, because for θ ≤ 5.09° the Γ-valley valence
bands lie below the 6th valence band, so Table I contains no Γ parameters.  The
defaults above reuse the AB values: at ψ_Γ = 0 the AA (B38) and AB (C13)
potentials are the same expression, so this is the least arbitrary available
choice.  They are **not** fits.

**2. `E_Γ` is not tabulated anywhere.**  (C19) is
`h = ħ²∇²/(2m_Γ*) + w_Γ + E_Γ + V_Γ(r)`, in which `E_Γ` only slides the whole
Γ-valley manifold relative to the ±K valleys.  App. B 3 says as much — the
zeroth harmonic "must be included to compare the energies of the in-active Γ
valley and active ±K valleys" — but never gives it, because each valley is
fitted separately.  It is therefore a free parameter, and the default is
*calibrated*, not fitted, from the paper's own statement (Sec. III.C and the
Fig. 4(c) caption) that the first ultra-flat Γ band lies 30 meV below the ±K
valence band maximum.  `E0_relative_to_vbm` (in `bands.jl`) recomputes it,
`03_selftest.jl` Part E checks that the tabulated number is reproduced to 3×10⁻⁴
meV, and since the Hamiltonian is linear in `E_Γ` only the *position* depends on
it — never the band structure or the flat-band splitting.  (The recomputation is
only reproducible at the truncation the stored value was made with: the Γ flat
band itself still moves 0.37 meV between Nshell = 5 and 7, so at Nshell = 5 the
rule returns −639.76 rather than −640.13.  That is a cutoff effect, not a
geometry one — `a0` and θ, the two things the guard compares, do not depend on
the cutoff.)

For AA stacking the same rule would put the Γ valley *inside* the Fig. 4(a)/(b)
window, contradicting the paper's statement that it lies below the 6th valence
band, so `AA_GAMMA` is placed 160 meV below the ±K valence band maximum
(−167 / −160 meV for the FH / FH+SH fits).  **That offset is a placeholder**,
flagged here and in `parameters.jl`.

---

## How correctness is judged

No comparison against Fig. 4 is used.  The criteria are the paper's own
identities and the appendices' own analytic results, and they are all in
`scripts/03_selftest.jl`:

**A. Operator-level identities** (`symmetry.jl`, residual ‖LHS−RHS‖ at a generic
k unless noted; all pass):

| identity | equation | models |
|---|---|---|
| Hermiticity | (B13) | all |
| `H(k + b_{M,i}) = H(k)` / Q-set closure | (B23) | ±K |
| C₃ covariance `P₃ H(C₃k) P₃† = H(k)` | (B10) | all |
| C₃ reps at Γ_M ∈ {1, ω, ω²} | (B10) | all |
| effective inversion `P H(−k) P† = H(k)` | (B18) / (B39) | AA±K, AA-Γ, AB-Γ (ε = 0) |
| mirror-T, unitary part | (B10)/(B14)/(C15) | AA±K, AA-Γ (C₂ᵧ𝒯); AB-Γ (C₂ₓ𝒯) |
| valley TR `P H_{−η}(k) P† = conj H_η(−k)` | (B10) | ±K |
| AB effective TR `P conj H(k) P† = H(−k)` | (C5) / (C14) | AB±K, AB-Γ (ε = 0) |
| `H_{+1}(k) = H_{−1}(k)` identically | (C5)+(C8) | AB±K |
| `spec H_{+1}(k) = spec H_{−1}(k)` | (B18) | AA±K |
| tunneling is strictly *inter*layer | (B13) | ±K |
| H exactly block diagonal in the layer (w = 0) | Table II | AB±K |
| `spec H_t(k) = spec H_b(C₂ₓk)` | (B14) | AB±K |
| exact 2-fold spin degeneracy | every term ∝ s₀ | Γ |

Two of these are *predicted failures* and are reported as such rather than
asserted: the effective inversion for AB±K (it is not a symmetry — (B18)), and
the AB effective 𝒦 for AA.  In addition, whenever `ε ≠ 0`, every identity whose
operator exchanges the layers **must** break — (B19)/(B40)/(C15) make the
displacement field `(−1)^l ε/2` flip sign under the exchange.  Part C drives all
of them at ε = 10 meV and confirms they do (1.5×10² – 1.1×10³).  This is not
decoration: at ε = 0 the Γ-valley models are layer-exchange symmetric *on their
own*, so a layer-flip error in a symmetry permutation passes either way, and
ε ≠ 0 is the only lever that exposes it.

**B. Main source ↔ cross source.**  `kmodels.jl` assembles `h_QQ'(k)` from the
analytic Fourier components of the *real-space* potentials; `transcribe.jl`
re-derives the same matrix by transcribing the printed *momentum-space*
equations (B27)/(C8)/(C18) verbatim, scanning every pair of basis states with
floating-point delta functions and never touching the integer-key or
Fourier-table machinery.  The two agree to `0.00e+00` for all five parameter
sets at both valleys.  (App. B 3 prints no momentum-space Γ form for AA, so that
model has no second source; it is covered by Part A instead.)

**C. Quantitative cross-checks.**  AA ±K: η = +1 and η = −1 have identical
spectra point by point.  AB ±K: with w = 0 the layer degeneracy is exact exactly
where C₂ₓk stays in the same C₃ orbit — along Γ_M–M_M (where C₂ₓ = C₃²), not
along Γ_M–K_M — which is the paper's Fig. 4(c) pattern of fourfold lines on the
60° ray and twofold at K_M.  Shell convergence: the top 6 ±K valence bands move
by ≤ 3×10⁻⁵ meV (AA FH 9×10⁻¹², FH+SH 6×10⁻⁸, AB 1×10⁻¹⁰) between Nshell = 4, 5,
6, with the basis dimension printed alongside so that a silent zero is visible.
`ħω` from (C25) reproduces the paper's 56.76 meV to 0.0014 meV.

**D. The Γ-valley flat bands vs App. C 3.**  See the next section.

**E. The parameter records and the validity guard** (Part G of the self-test,
which runs after A–F).  It freezes every fitted number and the whole band
structure — three Σ|E| checksums along Γ_M–K_M–M_M–Γ_M, the ±K valence band
maximum, ħω, the anharmonic shift — and asserts they are *bit-identical* to the
values from before the `ParamSet` refactor, because provenance metadata has to be
physically inert.  It then exercises the registry, each of the three mismatch
failures, both escape hatches, and the re-stamping of the `E_Γ` calibration.

The figures are a product, not a criterion.

---

## The Γ-valley flat bands, and the paper's 56.74 meV

App. C 3 reduces the AB Γ-valley model to an array of decoupled harmonic
oscillators: `m_Γ* = 10 m_e` makes the oscillator length `l = (ħ²/2m_Γ*k) / ħω`
≈ 3.66 Å, tiny inside the 52 Å moiré cell, so the top valence states are the
`(n_x, n_y)` levels of `−V_Γ(r)` expanded about the AA site, with
`ħω = ħ√(3V_Γ|b_{M,1}|²/m_Γ*) = 56.76 meV` (C25).  (C24) then gives

* `(n_x,n_y) = (0,0)`: **one** spinless state, spinless C₃ = 1  → 2-fold with spin,
* `(1,0),(0,1)`: **two** spinless states, C₃ = e^{±i2π/3}  → 4-fold with spin.

The numerical diagonalisation reproduces exactly that (and the paper's Fig. 4(c)
caption): at Γ_M the top level is single with C₃ = 1, the next is a doublet with
C₃ = e^{±i2π/3}, the two manifolds sit at −30.00 and −84.90 meV, and their
bandwidths over the whole MBZ are 0.003 and 0.011 meV — "extremely flat" in the
paper's words, because the oscillator states are so localized.

**The splitting is 54.84 meV, not the 56.74 meV the paper quotes.**  That
difference is not an error in either calculation; it is the quartic
anharmonicity of the lattice potential, which (C23) drops.  Expanding (C20)
exactly,

```
−V_Γ(r) = −6V_Γ + (3/2)V_Γ|b_{M,1}|² r² − (3/32)V_Γ|b_{M,1}|⁴ r⁴ + O(r⁶),
```

so the perturbation about the harmonic Hamiltonian (C23) is
`W = −(3/32)V_Γ|b|⁴r⁴`.  With the 2D oscillator moments `⟨r⁴⟩_(0,0) = 2l⁴` and
`⟨r⁴⟩_avg[(1,0),(0,1)] = 6l⁴`, first-order perturbation theory gives

```
ΔE₁ − ΔE₀ = −(3/8) V_Γ |b|⁴ l⁴ = −(1/4) ħ²|b_{M,1}|²/(2m_Γ*) = −1.865 meV,
```

independent of `V_Γ`.  So the exact lattice splitting should be
`56.761 − 1.865 = 54.897 meV`; the numerical diagonalisation gives 54.863 meV at
Nshell = 7 and 54.839 meV at Nshell = 9, i.e. the analytic prediction is
reproduced to 0.034 meV (`03_selftest.jl` Part E).  (`02_bands_AB.jl` prints
54.895 meV: it takes max − min over the whole path rather than at Γ_M, so it also
picks up the 0.011 meV dispersion.)  The paper's own 56.74 meV
agrees instead with (C25)'s harmonic `ħω` to 0.02 meV, consistent with App. C 3
describing the whole flat-band analysis as an approximation ("can be
approximated by an array of decoupled harmonic oscillators").  `anharmonic_shift`
in `parameters.jl` holds the formula, and its docstring holds this argument.

---

## Notes, conventions and caveats

**The Γ-valley label carries no layer information.**  In the ±K valley the
plane-wave label `Q = G_M + η(−1)^l q₁` *encodes* the layer, so `Q → −Q` flips
it automatically and (B18) needs nothing extra.  The Γ-valley label set is just
`{G_M}`, so wherever the paper's operation exchanges layers the permutation must
do it explicitly: (B39) `𝓘ψ†_{r,l}𝓘⁻¹ = ψ†_{−r,l̄}` for the AA effective
inversion, (C14) `𝒦ψ†_{r,l}𝒦⁻¹ = ψ†_{r,l̄}` for the AB effective TR, and (B14)
for C₂ᵧ𝒯.  `basis_perm(::GammaBasis, ...; layerswap = true)` is that switch.

**Mirror × T acts on the momentum label as the *other* reflection.**  The
composite `C₂ᵧ𝒯` is `f ↦ conj(f(C₂ᵧr))`; conjugating `e^{−ip·r}` reflects the
label a second time, so `C₂ᵧ𝒯` maps the label by C₂ₓ and `k → C₂ₓk`, while AB's
`C₂ₓ𝒯` maps it by C₂ᵧ.  This is forced, not chosen: only that reflection takes
`q₁ → −q₁`, which lies in the *other* family of (B23) and is therefore exactly
the layer flip (B10)/(B14) requires.  (`c2x_key_g` / `c2y_key_g` in
`symmetry.jl`.)

**App. B and App. C disagree about which layer is `l = 0`.**  (B17) uses
`l = t, b ↔ (−1)^l = +1, −1`, while (B36)/(B38) fix `(−1)^b = +1`, `(−1)^t = −1`.
This is only a naming choice, equivalent to ε → −ε; at ε = 0 — where every band
structure here is computed — it is irrelevant, and at ψ_Γ = 0 even the
potential is unaffected.  `layer_sign_kvalley` and `layer_sign_gamma` in
`parameters.jl` keep each appendix's own printed convention.

**Truncation.**  Both sets are cut in a way their own point group leaves
invariant, so the cut never breaks a symmetry.  The Γ-valley triangular set
`{G_M}` uses the hex shell `|m1|, |m2|, |m1+m2| ≤ Nshell` (the idea borrowed from
the older example script), which is C₆-symmetric by construction.  The ±K set
`{G_M ± q₁}` is *centred on ±q₁*, not on G_M = 0, so the same hex shell would
**not** be C₃-invariant — C₃ maps `G + q₁` to `(C₃G − b_{M,2}) + q₁`, and that
label can fall outside the shell.  It is therefore truncated by the disk
`|Q| ≤ Rcut` instead, which is invariant under C₃ and inversion by construction;
`Rcut` is inherited from the hex shell so that `Nshell` still controls the size.
The ±K models converge by Nshell = 5–6 (see Part D of the self-test).  The
Γ-valley models need **Nshell ≥ 7**: `m_Γ* = 10 m_e` makes the oscillator states
far more localized than the ±K states, and the splitting converges slowly from
above — 62.23 (N = 4), 56.53 (5), 55.09 (6), 54.86 (7), 54.84 (8), 54.84 (9)
meV.  `scripts/02_bands_AB.jl` therefore uses Nshell = 7.

**The *low* Γ-valley bands are not converged, and are not shown.**  With the
valence-band sign the Γ-valley spectrum is unbounded below, and the
antibonding/large-|G| states below roughly −300 meV keep moving as Nshell grows;
their width increases with Nshell rather than settling.  Only the *top* of the
Γ manifold — the flat bands discussed above — is a converged, basis-independent
statement.  The plotted window (0 … −100 meV) and the flat-band claims are
entirely inside that converged region; no other Γ-valley quantity is quoted
anywhere in this package.

**No OCR error in (B27) or (C18).**  An earlier plan for this work asserted two:
that (B27) has the `±g_i` deltas of (B21) exchanged, and that (C18)'s
`δ_{G,G'∓g_i}` differs in sign from (C13).  Both were checked against the
printed text and are wrong — (B27) and (C18) are each *consistent* with the
real-space equations through the matrix-element rule derived above, which is
also why `kmodels.jl` and `transcribe.jl` agree exactly.  The other OCR
observations in that plan are accurate: Table II's columns are misaligned in the
extracted text, and several equation numbers (B12, C22) are swallowed by the
extraction.

---

## Reference

Jia, Zhang, Yu, Liu, Song, Wang, Yao, *Moiré fractional Chern insulators. I.
Origin of the insulating states and quantum anomalous Hall effect*,
Phys. Rev. B **109**, 205121 (2024).  Equation numbers in parentheses, e.g.
(B27), refer to its appendices B and C.
