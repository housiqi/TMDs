# Hubbard Model Physics in Transition Metal Dichalcogenide Moir´e Bands

Fengcheng Wu,<sup>1</sup> Timothy Lovorn,<sup>2</sup> Emanuel Tutuc,<sup>3</sup> and A. H. MacDonald<sup>2</sup> 

<sup>1</sup>Materials Science Division, Argonne National Laboratory, Argonne, Illinois 60439, USA 

<sup>2</sup>Department of Physics, University of Texas at Austin, Austin, Texas 78712, USA <sup>3</sup>Department of Electrical and Computer Engineering, Microelectronics Research Center, The University of Texas at Austin, Austin, Texas 78758, USA 

(Received 8 April 2018; published 10 July 2018) 

Flexible long period moir´e superlattices form in two-dimensional van der Waals crystals containing layers that differ slightly in lattice constant or orientation. In this Letter we show theoretically that isolated flat moir´e bands described by generalized triangular lattice Hubbard models are present in twisted transition metal dichalcogenide heterobilayers. The hopping and interaction strength parameters of the Hubbard model can be tuned by varying the twist angle and the three-dimensional dielectric environment. When the flat moir´e bands are partially filled, candidate many-body ground states at some special filling factors include spin-liquid states, quantum anomalous Hall insulators, and chiral d-wave superconductors. 

DOI: 10.1103/PhysRevLett.121.026402 

Introduction.—Long-period superlattices form when two-dimensional crystals are overlaid with a small difference in lattice constant or orientation. When the twodimensional crystals are semiconductors or semimetals, their low-energy electronic degrees of freedom $( \mathrm { d . o . f . } )$ can [1] be accurately described using continuum models in which commensurability between the moir´e pattern and the atomic lattice plays no role. Because the continuum model Hamiltonians are periodic in space, their single-particle eigenstates satisfy Bloch’s theorem and form bands in momentum space, referred to as moir´e bands. The moir´e band Hamiltonian acts in a spinor space whose dimension is determined by the number of low-energy bands in the host two-dimensional crystal. In the moir´e band model of twisted bilayer graphene, for example, there are sixteen low-energy bands corresponding to the bilayer’s four triangular sublattices, and to spin and valley. Mott insulators and superconductivity have recently been discovered in the flat bands of twisted bilayer graphene [2–4]. In this Letter, we construct moir´e band Hamiltonians for holes in twisted heterobilayers formed from semiconducting transition metal dichalcogenides (TMDs), which have only two low-energy valence bands when the chemical potential is within the topmost valence bands as illustrated in Fig. 1(a), and therefore map to single-band Hubbard models. We show theoretically that isolated flat Bloch bands described by generalized triangular lattice Hubbard models are present in TMD heterobilayers, and that spin-liquid states are likely to occur when these bands are close to half-filling. The moir´e bilayers provide a new solidstate platform to simulate the Hubbard model, one in which model parameters such as band width, interaction strength, and band filling are widely tunable. 

This two-dimensional platform can be studied at accessible temperatures using a variety of experimental techniques, for example, transport and scanning tunneling microscopy. 

![](images/e11276b41b115be9b57f32917c37d7bfcda61c7383e69dddc0516a131acc144e.jpg)



(c)


![](images/11678033889d17841370f2f8c354bea3d424b2ef5d6471021d21ed7760485e77.jpg)



(d)


![](images/1eac7be272e528db5aa9874665ccd1795e8821db615bb34ff3880126ba5c9803.jpg)


![](images/0d6b2597c6282f7d9b7a15915455f0ba047fee17ae1398ef68857c5328244672.jpg)



FIG. 1. (a) Schematic band structure of monolayer $\mathrm { W S e } _ { 2 }$ with a large (small) spin splitting at valence (conduction) band extrema located at $\mathrm { t h e } \pm K$ valleys. (b) AA stacked $\mathrm { W } X _ { 2 } / \mathrm { M o } X _ { 2 }$ bilayers with an additional in-plane displacement d, and a twist angle θ. a and $\pmb { a } _ { 2 }$ are primitive translation vectors of $\mathsf { W } X _ { 2 } .$ . (c) Dependence of the ${ \mathrm { W S e } } _ { 2 }$ valence band maximum $\Delta ( d )$ on displacement d in AA stacked $\mathrm { W S e } _ { 2 } / \mathrm { M o S e } _ { 2 }$ with zero twist angle. $\Delta ( d )$ has triangular lattice periodicity and one maximum per triangular lattice unit cell. (d) When a moir´e pattern is formed the band maximum variation is magnified from the atomic scale to the moir´e pattern scale. The color scales in (c) and (d) are identical and the orange dashed lines in (d) are near-neighbor links that connect $\Delta ( r )$ maxima.


To simulate the two-dimensional Hubbard model, it is necessary to identify a moir´e system in which it is possible to isolate a single band with a twofold degeneracy. Because two-dimensional group-VI TMDs such as $\mathbf { M o S } _ { 2 }$ have band extrema at two inequivalent Brillouin-zone corners [5], conduction band states occur in groups of four with valley degeneracy and a small spin splitting due to spin-orbit interactions. We therefore focus on the isolation of moir´e bands formed from orbitals in the valence band, which have a very large spin splitting [Fig. 1(a)]. Valley d.o.f. then faithfully plays the role of spin in the Hubbard model. To avoid any additional degeneracy due to the presence of two layers, we choose to study heterobilayers. We consider common-chalcogen TMDs $\mathrm { W } X _ { 2 } / \mathrm { M o } X _ { 2 } \quad ( X = \mathrm { S } .$ Se), which have very similar lattice constants and can have long period moir´e patterns. Below we focus on the influence of the moir´e pattern on states near the maximum of the ${ \mathrm { W S e } } _ { 2 }$ valence bands, which lie inside the ${ \bf M o S e } _ { 2 }$ gaps and are only weakly coupled to states in ${ \bf M o S e } _ { 2 }$ due to band offsets [6]. We note that similar physics can also be realized in TMD bilayers with different chalcogen atoms [7,8], for example, $\mathrm { W S e } _ { 2 } / \mathrm { M o S } _ { 2 } .$ , which is studied in detail in the Supplemental Material [9]. 

Moir´e potential.—To derive the valence band moir´e Hamiltonian from first principles we follow the approach outlined in Ref. [11], which in the present case requires an evaluation of the dependence of the $\mathrm { W } X _ { 2 }$ valence band maximum energy $\Delta$ on the relative displacement d between two layers with identical lattice constants and twist angle $\theta = 0$ . The ab initio calculation was performed using fully relativistic density-functional theory in the local-density approximation as implemented in QUANTUM ESPRESSO [12]. In Fig. 1(c), we plot numerical values of $\Delta ( d )$ for the AA stacked $\mathrm { W S e } _ { 2 } / \mathrm { M o S e } _ { 2 }$ bilayer illustrated in Fig. 1(b). In the twisted bilayer moir´e pattern $( \theta \neq 0 )$ , the local value of d changes slowly over the moir´e period $( a _ { M } )$ and the valence band maximum, which serves as a spin-independent external potential, follows the variation of d, and varies periodically in space. Because we are interested only in moir´e periods greatly in excess of the host material lattice constant $( a _ { 0 } )$ , an effective mass approximation can be used for the band dispersion of the host material. We choose $m ^ { * } \sim 0 . 3 5 m _ { 0 }$ for $\mathrm { W S e } _ { 2 } ,$ where $m _ { 0 }$ is the free electron mass. Combining these considerations we obtain the following moir´e band Hamiltonian: 

$$
\begin{array}{c} \mathcal {H} = - \frac {\hbar^ {2} Q ^ {2}}{2 m ^ {*}} + \Delta (\boldsymbol {r}), \\ \Delta (\boldsymbol {r}) = \sum_ {\boldsymbol {b}} ^ {\prime} V (\boldsymbol {b}) \exp [ i \boldsymbol {b} \cdot \boldsymbol {r} ], \end{array}\tag{1}
$$

where $- \hbar ^ { 2 } Q ^ { 2 } / ( 2 m ^ { * } )$ and $\Delta ( r )$ are the moir´e band kinetic and potential energies. The potential $\Delta ( r )$ shown in Fig. 1(d) can be accurately approximated by a Fourier expansion that includes only the six moir´e reciprocal lattice vectors in the first shell. Because the potential is real and each TMD monolayer has threefold-rotational symmetry, we require $V ( \pmb { b } ) = V ^ { * } ( - \pmb { b } )$ and $V ( \hat { \mathcal { R } } _ { 2 \pi / 3 } \pmb { b } ) = V ( \pmb { b } )$ Therefore, all six $V ( b )$ are fixed by $V ( pmb { b } _ { 1 } ) = V \exp ( i \psi )$ where $\pmb { b } _ { 1 } = 4 \pi / ( \sqrt { 3 } a _ { M } ) \hat { x }$ . Fitting to the ab initio potential energy, we find that $( V , \psi )$ is $( 6 . 6 \mathrm { m e V } , - 9 4 ^ { \circ } )$ for ${ \mathrm { W S e } } _ { 2 }$ on $\mathbf { M o S e } _ { 2 }$ in AA stacking. The fitting procedure has been described in detail in Refs. [13,14]. Because the coupling between the two layers can be modified by external vertical electric field [15] and by pressure [16], the moir´e potential is experimentally tunable. The unprecedented advantage of van der Waals heterobilayers is that the moir´e potential period can be tuned simply by changing the twist angle: $a _ { M } \approx a _ { 0 } / \theta .$ . In the case of $\mathrm { W S e } _ { 2 } / \mathrm { M o S e } _ { 2 } , a _ { M }$ is about 19 nm at $1 ^ { \circ }$ twist angle. We note that collective excitations, for example, excitons, experience a similar moir´e potential energy whose influence has been studied theoretically in Refs. [13,14,17]. 

Hubbard model.—The length scales relevant to moir´e Hubbard band formation are the moir´e period $a _ { M }$ and the spatial extend $a _ { W }$ of the Wannier functions associated with the highest-energy moir´e band, which is localized around the triangular lattice of moir´e potential maximum positions. Near its maximum the moir´e potential can be approximated by a harmonic oscillator potential: $- \beta V ( \delta r / a _ { M } ) ^ { 2 } / 2$ , where $\beta =$ $1 6 \pi ^ { 2 } \cos ( \psi + 1 2 0 ^ { \circ } )$ for the potential shown in Fig. 1(d). Within this approximation, $a _ { W } \approx [ \hbar ^ { 2 } / ( \beta m ^ { * } V ) ] ^ { 1 / 4 } \sqrt { a _ { M } } .$ Because $a _ { W } / a _ { M }$ scales as $1 / \sqrt { a _ { M } }$ , we can anticipate that the highest energy moir´e band flattens with a decrease in the twist angle. 

In Fig. 2(a) we plot the moir´e bands of ${ \mathrm { W S e } } _ { 2 }$ on ${ \bf M o S e } _ { 2 }$ obtained by diagonalizing the moir´e Bloch Hamiltonian $\mathcal { H } ( k )$ in a plane wave representation: 

$$
\langle \boldsymbol {k} + \boldsymbol {g} ^ {\prime} | \mathcal {H} | \boldsymbol {k} + \boldsymbol {g} \rangle = - \delta_ {\boldsymbol {g} ^ {\prime}, \boldsymbol {g}} \frac {\hbar^ {2} | \boldsymbol {k} + \boldsymbol {g} | ^ {2}}{2 m ^ {*}} + V (\boldsymbol {g} ^ {\prime} - \boldsymbol {g}),\tag{2}
$$

where $\pmb { g }$ and $\pmb { g } ^ { \prime }$ are moir´e reciprocal lattice vectors. The highest valence moir´e band at $\bar { \theta } = 2 . 0 ^ { \circ }$ is separated from other bands by an energy gap and has a narrow bandwidth $\left( \sim 1 1 \ \mathrm { m e V } \right)$ . This isolated flatband can be described by a tight-binding model on a triangular lattice: 

$$
H _ {0} = \sum_ {\tau = \uparrow , \downarrow} \sum_ {\boldsymbol {R}, \boldsymbol {R} ^ {\prime}} t (\boldsymbol {R} ^ {\prime} - \boldsymbol {R}) c _ {\boldsymbol {R} \tau} ^ {\dagger} c _ {\boldsymbol {R} ^ {\prime} \tau},\tag{3}
$$

where R represents the triangular lattice formed by the moir´e potential maximum positions, and τ is a valley index. In Fig. 2(b), we show the density of states (DOS) of the single-particle moir´e bands as a function of hole density, which is strongly enhanced by the moir´e potential, and has sharp peaks at moir´e band saddle points. The flatband energy dispersion can be accurately fit by including hopping up to the third nearest neighbor. Figure 2(d) shows the hopping parameters $t _ { n }$ as a function of moir´e period $a _ { M }$ , where $t _ { n }$ connects the nth nearest neighbors. The hopping parameters are real, $\left| t _ { 1 } \right|$ is dominant over $\left| t _ { 2 , 3 } \right|$ and all three hopping parameters decrease exponentially with increasing $a _ { M }$ 

![](images/aa331ca3fec63b600405c731f1b6ac98b46531060241923c085ebd6438a9e884.jpg)


![](images/636bad5ccae49c087952ffaeca7d9d263d47cc423093cfd1e790e04395066885.jpg)



(c)


![](images/57b3724df45b257c380b0c9d674403b42e1fb95c2cdf11e3f999a80a745a337d.jpg)


![](images/c8c6787f2ca232112edd156a766eedfc4c39a9dc381e8ad3edb51081f496d14a.jpg)



FIG. 2. (a) Moir´e bands at twist angle $\theta = 2 . 0 ^ { \circ } .$ . The red dashed line is a tight-binding-model fit to the highest valence band that includes hopping up to the third nearest neighbor. (b) Density of states as a function of the hole filling factor $n _ { h } / n _ { M }$ (bottom) and the hole density $n _ { h }$ (top). (c) The Wannier function $w ( r )$ associated with the highest-energy band in (a). w r is centered on one of the moir´e potential maxima positions. The orange lines are the links of the moir´e triangular lattice. (d) Hopping parameters $t _ { n }$ vs neighbor number n as a function of moir´e period $a _ { M }$ (bottom) and twist angle θ (top).


Figure 2(c) plots the Wannier wave function $w ( r )$ constructed from the isolated band’s Bloch states. The spatial extent $a _ { W }$ of this localized wave function increases with moir´e pattern period, in agreement with the estimate above, but its ratio to $a _ { M }$ decreases. Correspondingly the on-site Coulomb repulsion energy $U _ { 0 } \sim e ^ { 2 } / ( \epsilon a _ { W } )$ decreases slowly as the moir´e period increases. It follows that the ratio of $U _ { 0 }$ to the bandwidth increases very quickly with $a _ { M }$ , and that electronic states formed when the moir´e band is partially occupied by electrons become strongly correlated. The effective dielectric constant ϵ in the bilayer is sensitive to the three-dimensional dielectric environment out to vertical distances ${ \sim } a _ { M }$ from the bilayer, allowing the strength of correlations at a given orientation angle to be adjusted over a wide range. To simulate a Hubbard model with short-range repulsion, we assume that a metallic screening layer is close to the TMD bilayer, but separated from it by a dielectric. Such a metallic layer, formed by graphene, for example, could also act as a gate that controls the filling factor of the moir´e band. In a simple imagecharge approximation, the electron-electron interaction potential is ${ \tilde { U } } ( { \boldsymbol { r } } ) = ( e ^ { 2 } / \epsilon ) [ r ^ { - 1 } - ( r ^ { 2 } + D ^ { 2 } ) ^ { - 1 / 2 } ]$ where $D / 2$ is the vertical distance between the metallic layer and the TMD bilayer. When $\tilde { U } ( \boldsymbol { r } )$ is projected onto the isolated band Wannier states, and the negligible overlap between Wannier orbitals centered on different sites is noted, the interaction Hamiltonian reduces to the generalized Hubbard form 

$$
H _ {1} = \frac {1}{2} \sum_ {\tau , \tau^ {\prime}} \sum_ {\boldsymbol {R}, \boldsymbol {R} ^ {\prime}} U (\boldsymbol {R} ^ {\prime} - \boldsymbol {R}) c _ {\boldsymbol {R}, \tau} ^ {\dagger} c _ {\boldsymbol {R} ^ {\prime}, \tau^ {\prime}} ^ {\dagger} c _ {\boldsymbol {R} ^ {\prime}, \tau^ {\prime}} c _ {\boldsymbol {R}, \tau}.\tag{4}
$$

In Fig. 3 we plot values of repulsive interaction $U _ { 0 }$ (on site), $U _ { 1 }$ (nearest neighbor), and $U _ { 2 }$ (second nearest neighbor) as a function of moir´e period $a _ { M }$ . These results were calculated using $D = 3$ nm, and $U _ { 1 , 2 }$ are therefore strongly suppressed compared to $U _ { 0 }$ . As expected from the scaling analysis above, $U _ { 0 }$ decreases only slowly as $a _ { M }$ increases. 

Equation (3) combined with (4) describes a generalized Hubbard model on a triangular lattice. The isolated band is fully occupied when the TMD bilayer is charge neutral. By inducing hole carriers, the band becomes partially occupied. When the isolated band is completely emptied by hole doping, the carrier density is $\hat { n _ { M } } = \bar { 2 / ( \sqrt { 3 } \hat { a } _ { M } ^ { 2 } / 2 ) }$ , where the factor of 2 accounts for the Hubbard model spin degeneracy. We find that $n _ { M } = ( 0 . 6 4 \times 1 0 ^ { 1 2 } ) ( \theta [ ^ { \circ } ] ) ^ { 2 } ~ \mathrm { c m } ^ { - 2 }$ implying that the full range of band fillings is accessible by electrical gating for θ less than ${ \sim } 4 ^ { \circ }$ . In the following, we discuss possible moir´e band ground states at $1 / 2$ and $3 / 4$ hole doping. 

Halffilling.—When the isolated band is halffilled, there is one electron per moir´e unit cell. As illustrated in Fig. 3(a), $U _ { 0 } \gg t _ { 1 }$ is satisfied even for a relatively large dielectric constant employed to obtain these estimates. The strong onsite repulsion suppresses double occupation of moir´e lattice sites and gives rise to a Mott insulator ground state with only spin d.o.f. at low energies. In the large $U _ { 0 }$ limit, the Hubbard model can be mapped to the spin Heisenberg model: 

![](images/33cebc3f09e3c26bacf71c43202a23d8a270244732ac1965b3ab438d1d3460f0.jpg)



1.0°


![](images/abd420a30c9d6f92ce2d54fe6dabb4eb106b5f5677d923c4deda7b6cebd4796a.jpg)



FIG. 3. (a) Hubbard model repulsive interaction parameters $\epsilon U _ { n } .$ , and (b) spin exchange interactions $J _ { n }$ as a function of moir´e period $a _ { M }$ and twist angle $\theta .$ The insets of (a) and (b), respectively, show the ratios $U _ { 0 } / t _ { 1 }$ and $J _ { 2 } / J _ { 1 }$ . The semiconductor background dielectric constant ϵ was set to 10 for the evaluation of $U _ { n } .$ Smaller effective values of ϵ are applicable when the dielectric environment is engineered to maximize interaction strength. For example, ϵ is about 5 when hexagonal boron nitride is used as the dielectric layer [18].


$$
H _ {S} = \sum_ {\boldsymbol {R}, \boldsymbol {R} ^ {\prime}} ^ {\prime} J (\boldsymbol {R} ^ {\prime} - \boldsymbol {R}) \mathbf {S} _ {\boldsymbol {R}} \cdot \mathbf {S} _ {\boldsymbol {R} ^ {\prime}},\tag{5}
$$

where S is the $S = 1 / 2$ spin operator, J is a spin exchange coupling energy, and the prime on the sum indicates that each pair of sites is counted only once. Using $t / U$ perturbation theory [19] to calculate the exchange interactions up to the third nearest neighbors, we find that $J _ { 1 } = 4 ( t _ { 1 } ^ { 2 } / U _ { 0 } ) [ 1 - 7 ( t _ { 1 } / U _ { 0 } ) ^ { 2 } ] , J _ { 2 } = 4 t _ { 2 } ^ { 2 } / U _ { 0 } + 4 t _ { 1 } ^ { 4 } / U _ { 0 } ^ { 3 } ,$ and $J _ { 3 } = 4 t _ { 3 } ^ { 2 } / U _ { 0 } + 4 t _ { 1 } ^ { 4 } / U _ { 0 } ^ { 3 }$ . Here we have expanded to second order in $t _ { 2 , 3 }$ , but to fourth order in $t _ { 1 }$ because $\left| t _ { 1 } \right| \gg \left| t _ { 2 , 3 } \right|$ . The numerical values of $J _ { n }$ are plotted in Fig. 3(b). 

The properties of triangular-lattice Heisenberg models have been thoroughly investigated in previous work. When only nearest-neighbor coupling $J _ { 1 }$ is nonzero, the ground state has three-sublattice $1 2 0 ^ { \circ }$ long range antiferromagnetic order. This antiferromagnetic state becomes unstable when the second nearest-neighbor coupling $J _ { 2 }$ exceeds a critical value. For the quantum spin-1=2 Heisenberg model on triangular lattice, a spin liquid phase has been found in the parameter region $0 . 0 \bar { 6 } \lesssim J _ { 2 } / J _ { 1 } \lesssim 0 . 1 7$ [20,21]. As shown in Fig. 3(b), $J _ { 2 } / J _ { 1 }$ exceeds 0.06 when the twist angle θ is larger than 3.0°, which makes the spin liquid state likely to occur. In our case, $J _ { 3 }$ is also nonzero but its small magnitude seems unlikely to significantly alter earlier estimates of phase boundaries. The arrival of moir´e band strong correlation physics motivates new studies of Heisenberg models with exchange coupling to further neighbors. 

3/4 filling.—At 3=4 hole doping, flatbands with only nearest-neighbor hopping have a van Hove singularity and the corresponding Fermi surfaces are perfectly nested. As shown in Figs. 2(b) and 4(a), these features are largely retained in realistic models because remote neighbor hopping is weak. The nesting vectors are $\pmb { b } / 2 ,$ , where b is a first-shell moir´e reciprocal lattice vector. One candidate ground state has the four-sublattice tetrahedral antiferromagnetic order [22] illustrated in Fig. 4. This magnetic order fully gaps the nested Fermi surface, and gives rise to a quantum anomalous 

![](images/76ed861d3a208aa03239caf2dd6985e1c8b749d2a59803768408cb1659bb2d48.jpg)



(b)


![](images/df99cf9483f7dde6ad06d4586c623cc1ee56cccd151a0c2b70c96d98a68e000b.jpg)



(c)


![](images/bb42af55b184c4ffab87aeda5701e7d53896621eccfd61e291c7aada5ce32b3f.jpg)



FIG. 4. (a) Energy contours in the moir´e Brillouin zone for the highest energy band in Fig. 2(a). (b) Magnetic order on a triangular lattice with four-sublattice tetrahedral antiferromagnetic order. The magnetic moment directions on the four sublattices are specified by the corresponding arrows in (c).


Hall insulator [22] with quantized Hall conductivity of $e ^ { 2 } / h$ at $T = 0$ . The tetrahedral order is noncoplanar, and results in a scalar spin chirality: $\chi = S _ { i } \cdot ( S _ { j } \times S _ { k } )$ . Thermal fluctuations at finite temperature will destroy long-range magnetic order in two dimension. However, the chirality χ is an Ising order parameter which can persist even at finite temperatures and support an anomalous Hall effect. Another candidate state has four-sublattice collinear antiferromagnetic order with site-dependent spin moments [23]. This state has gapless charge excitations at the Fermi energy in one spin component only, and therefore is a half metal that supports spin currents. In close competition with these magnetic states, there is also an instability towards chiral d-wave superconductivity from repulsive interactions [24,25]. In the renormalization group analysis, d-wave superconductivity has been found to be the leading weak coupling instability at $3 / 4$ filling [25]. 

Discussion.—For twist angles smaller than around $3 . 5 ^ { \circ } .$ the highest energy $\mathrm { W S e } _ { 2 }$ valence moir´e band provides a realization of the triangular lattice Hubbard models. For the special case of half-filling the system provides a realization of quantum spin models on triangular lattices. Although the triangular lattice is frustrated, the spin-model ground state is a relatively conventional antiferromagnet when only nearest-neighbor interactions are present. Our calculations demonstrate that twist angles can be turned to regimes in which spin-liquid states are expected. The estimated spininteraction energy scales are on the millielectron volt energy scale, making the low-temperature properties of these quantum spin systems accessible at dilution fridge temperature scales. The competition between strongly correlated states in the moir´e band Hubbard model can be tuned by the twist angle, the dielectric environment, and by strain that generates anisotropy for the triangular lattice. Furthermore, moir´e band Hubbard model realizations also allow strongly correlated electron systems to be studied in new ways. For example by examining how the carrier density depends on gate voltages it is possible to extract the Hubbard model chemical potential as a function of carrier density, and in this way to quantitatively extract, among other properties, the size of charge gaps expected at $1 / 2$ filling, and in some cases also at $3 / 4$ filling. 

One of the most interesting possibilities offered by TMD moir´e band systems is that of measuring spin transport characteristics in strongly correlated electron systems and comparing them directly with charge transport characteristics. Single layer TMD systems can be optically driven [26–28] into steady states with valley (and therefore spin) dependent chemical potentials, allowing them to be used as spin and charge reservoirs, and as spin-polarization detectors. These capabilities allow for measurements of coupled spin and charge transport in strongly correlated electron systems, a topic of great theoretical interest [29,30], and an important goal of cold-atom Hubbard model simulation efforts [31–33]. 

F. W. thanks I. Martin for valuable discussions. Work at Austin was supported by the Department of Energy, Office of Basic Energy Sciences under Contract No. DE-FG02- 02ER45958 and Award No. DE-SC0012670, and by the Welch Foundation under Grant No. TBF1473. Work at Argonne National Laboratory was supported by the Department of Energy, Office of Science, Materials Science and Engineering Division. The authors acknowledge HPC resources provided by the Texas Advanced Computing Center (TACC) at The University of Texas at Austin. 



[1] R. Bistritzer and A. H. MacDonald, Proc. Natl. Acad. Sci. U.S.A. 108, 12233 (2011). 





[2] K. Kim, A. DaSilva, S. Huang, B. Fallahazad, S. Larentis, T. Taniguchi, K. Watanabe, B. J. LeRoy, A. H. MacDonald, and E. Tutuc, Proc. Natl. Acad. Sci. U.S.A. 114, 3364 (2017). 





[3] Y. Cao, V. Fatemi, A. Demir, S. Fang, S. L. Tomarken, J. Y. Luo, J. D. Sanchez-Yamagishi, K. Watanabe, T. Taniguchi, E. Kaxiras, R. C. Ashoori, and P. Jarillo-Herrero, Nature (London) 556, 80 (2018). 





[4] Y. Cao, V. Fatemi, S. Fang, K. Watanabe, T. Taniguchi, E. Kaxiras, and P. Jarillo-Herrero, Nature (London) 556, 43 (2018). 





[5] D. Xiao, G.-B. Liu, W. Feng, X. Xu, and W. Yao, Phys. Rev. Lett. 108, 196802 (2012). 





[6] C. Zhang, C. Gong, Y. Nie, K.-A. Min, C. Liang, Y. J. Oh, H. Zhang, W. Wang, S. Hong, L. Colombo, R. M. Wallace, and K. Cho, 2D Mater. 4, 015026 (2017). 





[7] C. Zhang, C.-P. Chuu, X. Ren, M.-Y. Li, L.-J. Li, C. Jin, M.-Y. Chou, and C.-K. Shih, Sci. Adv. 3, e1601459 (2017). 





[8] Y. Pan, S. Fölsch, Y. Nie, D. Waters, Y.-C. Lin, B. Jariwala, K. Zhang, K. Cho, J. A. Robinson, and R. M. Feenstra, Nano Lett. 18, 1849 (2018). 





[9] See Supplemental Material at http://link.aps.org/ supplemental/10.1103/PhysRevLett.121.026402 for moir´e bands in $\mathrm { W S e } _ { 2 } / \mathrm { M o S } _ { 2 }$ bilayer and additional discussion on moir´e potential, which includes Ref. [10]. 





[10] T. Björkman, A. Gulans, A. V. Krasheninnikov, and R. M. Nieminen, Phys. Rev. Lett. 108, 235502 (2012). 





[11] J. Jung, A. Raoux, Z. Qiao, and A. H. MacDonald, Phys. Rev. B 89, 205414 (2014). 





[12] P. Giannozzi et al., J. Phys. Condens. Matter 21, 395502 (2009). 





[13] F. Wu, T. Lovorn, and A. H. MacDonald, Phys. Rev. Lett. 118, 147401 (2017). 





[14] F. Wu, T. Lovorn, and A. H. MacDonald, Phys. Rev. B 97, 035306 (2018). 





[15] H. C. P. Movva, T. Lovorn, B. Fallahazad, S. Larentis, K. Kim, T. Taniguchi, K. Watanabe, S. K. Banerjee, A. H. MacDonald, and E. Tutuc, Phys. Rev. Lett. 120, 107703 (2018). 





[16] M. Yankowitz, J. Jung, E. Laksono, N. Leconte, B. L. Chittari, K. Watanabe, T. Taniguchi, S. Adam, D. Graf, and C. R. Dean, Nature (London) 557, 404 (2018). 





[17] H. Yu, G.-B. Liu, J. Tang, X. Xu, and W. Yao, Sci. Adv. 3, e1701696 (2017). 





[18] C. R. Dean, A. F. Young, I. Meric, C. Lee, L. Wang, S. Sorgenfrei, K. Watanabe, T. Taniguchi, P. Kim, K. L. Shepard, and J. Hone, Nat. Nanotechnol. 5, 722 (2010). 





[19] A. H. MacDonald, S. M. Girvin, and D. Yoshioka, Phys. Rev. B 37, 9753 (1988). 





[20] Z. Zhu and S. R. White, Phys. Rev. B 92, 041105 (2015). 





[21] W.-J. Hu, S.-S. Gong, W. Zhu, and D. N. Sheng, Phys. Rev. B 92, 140403 (2015). 





[22] I. Martin and C. D. Batista, Phys. Rev. Lett. 101, 156402 (2008). 





[23] R. Nandkishore, G.-W. Chern, and A. V. Chubukov, Phys. Rev. Lett. 108, 227204 (2012). 





[24] R. Nandkishore, L. S. Levitov, and A. V. Chubukov, Nat. Phys. 8, 158 (2012). 





[25] R. Nandkishore, R. Thomale, and A. V. Chubukov, Phys. Rev. B 89, 144501 (2014). 





[26] K. F. Mak, K. He, J. Shan, and T. F. Heinz, Nat. Nanotechnol. 7, 494 (2012). 





[27] K. F. Mak, K. L. McGill, J. Park, and P. L. McEuen, Science 344, 1489 (2014). 





[28] K. Hao, G. Moody, F. Wu, C. K. Dass, L. Xu, C.-H. Chen, L. Sun, M.-Y. Li, L.-J. Li, A. H. MacDonald, and X. Li, Nat. Phys. 12, 677 (2016). 





[29] H. Kim and D. A. Huse, Phys. Rev. A 86, 053607 (2012). 





[30] C. Karrasch, D. M. Kennes, and J. E. Moore, Phys. Rev. B 90, 155104 (2014). 





[31] T. Esslinger, Annu. Rev. Condens. Matter Phys. 1, 129 (2010). 





[32] R. A. Hart, P. M. Duarte, T.-L. Yang, X. Liu, T. Paiva, E. Khatami, R. T. Scalettar, N. Trivedi, D. A. Huse, and R. G. Hulet, Nature (London) 519, 211 (2015). 





[33] M. A. Nichols, L. W. Cheuk, M. Okan, T. R. Hartke, E. Mendez, T. Senthil, E. Khatami, H. Zhang, and M. W. Zwierlein, arXiv:1802.10018. 

