"""
    MoireTMD

Continuum models of twisted bilayer TMDs, following PRB **109**, 205121 (2024)
(Jia et al., *Moire fractional Chern insulators...*).

Four single-particle models are implemented, exactly as given in the paper:

| model      | stacking | valley | equations                          |
|------------|----------|--------|------------------------------------|
| `:AA_K`    | AA       | +-K    | (B15), (B17), (B20), (B21) -> (B27)|
| `:AB_K`    | AB       | +-K    | (C4) -> (C8),  Eq. (5)             |
| `:AA_Gamma`| AA       | Gamma  | (B28)-(B38)                        |
| `:AB_Gamma`| AB       | Gamma  | (C11)-(C18),   Eq. (7)             |

Every model is diagonalised in the plane-wave basis of (B22)/(C16)
"""
module MoireTMD

using LinearAlgebra
using Printf

include("geometry.jl")
include("parameters.jl")
include("potentials.jl")
include("kmodels.jl")
include("transcribe.jl")
include("symmetry.jl")
include("bands.jl")

export MoireGeometry, moire_geometry, check_geometry, lattice_coeffs, is_lattice_vec,
       KValleyBasis, kvalley_basis, GammaValleyBasis, gammavalley_basis,
       GammaBasis, gammabasis, KPath, bz_path,
       KParams, GammaParams,
       AA_K_FH, AA_K_FHSH, AB_K, AB_GAMMA, AA_GAMMA,
       A0_MOTE2, THETA_DEG, HBAR2_2ME,
       kinetic_prefactor, hbar_omega, harmonic_stiffness, anharmonic_shift,
       layer_sign_kvalley, layer_sign_gamma,
       build_basis, hamiltonian, ishermitian_exact,
       ham_AA_K, ham_AB_K, ham_AA_Gamma, ham_AB_Gamma,
       transcribe_hamiltonian,
       BandData, diagonalize_path, diagonalize_point, top_bands, rigid_shift!, shifted,
       vbm, with_E0, E0_relative_to_vbm,
       CheckResult, symmetry_report, print_report, basis_perm, c3_rep, gidx,
       top_multiplets, c3_perm_spinless, key_family,
       c3_key, c3_key_inv, inv_key, c2x_key, pair_degeneracy, pair_gap,
       c3_key_g, c3_key_g_inv, inv_key_g, c2x_key_g, c2y_key_g,
       transcribe_AA_K, transcribe_AB_K, transcribe_AB_Gamma

end
