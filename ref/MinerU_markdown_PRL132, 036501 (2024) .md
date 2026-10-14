# Fractional Chern Insulator in Twisted Bilayer MoTe<sub>2</sub>

Chong Wang ,<sup>1</sup> Xiao-Wei Zhang,<sup>1</sup> Xiaoyu Liu,<sup>1</sup> Yuchi He ,<sup>2</sup> Xiaodong Xu,<sup>3,1</sup> Ying Ran,<sup>4</sup> Ting Cao,<sup>1,*</sup> and Di Xiao <sup>1,3,†</sup> 

<sup>1</sup>Department of Materials Science and Engineering, University of Washington, Seattle, Washington 98195, USA 

<sup>2</sup>Rudolf Peierls Centre for Theoretical Physics, Clarendon Laboratory, Parks Road, Oxford OX1 3PU, United Kingdom 

<sup>3</sup>Department of Physics, University of Washington, Seattle, Washington 98195, USA 

<sup>4</sup>Department of Physics, Boston College, Chestnut Hill, Massachusetts 02467, USA 

0 (Received 24 April 2023; accepted 30 November 2023; published 17 January 2024) 

A recent experiment has reported the first observation of a zero-field fractional Chern insulator (FCI) phase in twisted bilayer $\mathrm { M o T e } _ { 2 }$ moir´e superlattices [J. Cai et al., Signatures of fractional quantum anomalous Hall states in twisted $\begin{array} { r } { \mathbf { M o T e } _ { 2 } . } \end{array}$ , Nature (London) 622, 63 (2023).]. The experimental observation is at an unexpected large twist angle $3 . 7 ^ { \circ }$ and calls for a better understanding of the FCI in real materials. In this Letter, we perform large-scale density functional theory calculation for the twisted bilayer MoTe and find that lattice reconstruction is crucial for the appearance of an isolated flat Chern band. The existence of the FCI state at $\nu = - 2 / 3$ is confirmed by exact diagonalization. We establish phase diagrams with respect to the twist angle and electron interaction, which reveal an optimal twist angle of 3.5° for the observation of FCI. We further demonstrate that an external electric field can destroy the FCI state by changing band geometry and show evidence of the $\nu = - 3 / 5$ FCI state in this system. Our research highlights the importance of accurate single-particle band structure in the quest for strong correlated electronic states and provides insights into engineering fractional Chern insulator in moir´e superlattices. 

DOI: 10.1103/PhysRevLett.132.036501 

Introduction.—Fractional Chern insulators (FCIs), the analog of the fractional quantum Hall effect [1,2] in lattice systems, feature fractional excitations with anyonic statistics [3–6]. Because of its exotic nature and potential applications in topological quantum computing [7,8], the FCI has long been sought after in condensed matter experiments. To find the FCI, the widely accepted approach is to find a Chern band with significantly quenched kinetic energy [9–13]. This can be realized in two-dimensional (2D) moir´e superlattices, which have been shown to be a fruitful and tunable platform to explore electronic correlation [14,15]. Indeed, theoretical proposals of FCIs have been put forward in twisted bilayer graphene [16–21] and twisted bilayer transition metal dichalcogenide (TMD) [22–24]. Experimental evidences for FCIs has also been found in graphene-based superlattices [25,26], albeit at a finite magnetic field. Recently, the first observation of the FCI in the absence of magnetic field (i.e., the fractional quantum anomalous Hall effect) has been reported in twisted bilayer $\mathbf { M o T e } _ { 2 }$ at hole fillings $\nu = - 2 / 3$ and $\nu = - 3 / 5 \ [ 2 7 ]$ . The observation is at an unexpected large twist angle (∼3.7°), for which the Chern bands had been commonly believed to be too dispersive to stabilize the FCI. The experimental observation calls for a better understanding of single-particle band structure of twisted TMD bilayer as well as the emergence of the FCI in real materials. 

In this Letter, we perform large-scale density-functional theory (DFT) calculation for the twisted bilayer MoTe<sub>2</sub>. 

In contrast to previous theoretical studies, we find that the band structure features an isolated flat Chern band that favors FCIs. By exact diagonalization (ED) calculations, we confirm the existence of the FCI at $\nu = - 2 / 3$ . We also find that the ferromagnetism at $\nu = - 1 / 3$ is much weaker than that at $\nu = - 2 / 3$ , explaining the absence of the $\nu =$ −1=3 FCI in the experiment. We further investigate the fate of FCIs under an external electric field and find that FCIs become unstable at $E = 1 . 2 6 ~ \mathrm { m V / \mathring { A } }$ , which is consistent with experimental observation. The suppression of FCIs by external electric field is attributed to the deterioration of the flatness of band geometry. We establish phase diagrams with respect to the twist angle and electronic interaction, revealing an optimal twist angle of 3.5° for the observation of FCIs. Finally, we show evidences of the $\nu = - 3 / 5$ FCI state in this system. Our research highlights the importance of accurate single-particle band structure in the quest for strong correlated electronic states and provides insights into engineering FCIs in moir´e superlattices. 

Single-particle electronic structure.—Moir´e superlattices formed by twisting a bilayer introduce a long wavelength periodic structure characterized by the moir´e lattice constant $a _ { M } = a / \theta ;$ , where a is the lattice constant of the original 2D layer and θ is the twist angle. This large length scale makes it possible to model the low-energy electronic structure with a continuum model. The valence band edge of monolayer $\mathbf { M o T e } _ { 2 }$ is located at the corners of its Brillouin zone, i.e., K and K<sup>0</sup> points. The two points are separated by a large momentum such that K and K<sup>0</sup> valleys can be considered independently. Following the experiment [27], we consider R stacking twisted bilayer and the continuum model Hamiltonian for K valley reads [28,29] 


TABLE I. Parameters for the continuum model.


<table><tr><td></td><td>v (meV)</td><td>ψ (deg)</td><td>w (meV)</td></tr><tr><td>Local stacking approx. [28]</td><td>8.0</td><td>-89.6</td><td>-8.5</td></tr><tr><td>Large-scale DFT</td><td>20.8</td><td>+107.7</td><td>-23.8</td></tr></table>

$$
\mathcal {H} _ {K} = \left( \begin{array}{c c} H _ {b} & \Delta_ {T} (\boldsymbol {r}) \\ \Delta_ {T} ^ {\dagger} (\boldsymbol {r}) & H _ {t} \end{array} \right).\tag{1}
$$

Here, $H _ { b / t } = - \hbar ^ { 2 } ( k - K _ { b / t } ) ^ { 2 } / 2 m ^ { * } + \Delta _ { b / t } ( \pmb { r } )$ is the bottom (b) and top (t) layer Hamiltonian subjected to a moir´e potential $\begin{array} { r } { \bar { \Delta _ { b / t } } ( \pmb { r } ) = 2 v \sum _ { j = 1 , 3 , 5 } \cos \left( \pmb { G } _ { j } \cdot \pmb { r } \pmb { \bot } \psi \right) } \end{array}$ , where the bottom (top) layer corresponds to the positive (negative) sign. $K _ { b / t }$ is the K point for the bottom and top layer and $G _ { j }$ is the moir´e reciprocal lattice vectors defined by $G _ { j } =$ $( 4 \pi / \sqrt { 3 } a _ { M } ) \{ \cos [ \pi ( j - 1 ) / 3 ] , \sin [ \pi ( j - 1 ) / 3 ] \}$ . The interlayer tunneling is dictated by threefold rotational symmetry as $\Delta _ { T } ( r ) = w ( 1 + e ^ { - i G _ { 2 } \cdot r } + e ^ { - i G _ { 3 } \cdot r } )$ $m ^ { * } = 0 . 6 m _ { e }$ is the effective mass $( m _ { e }$ is the bare electron mass) and $( v , \psi , w )$ are the free parameters in the continuum model. The continuum Hamiltonian for the $K ^ { \prime }$ valley can be obtained by applying time-reversal symmetry to $\mathcal { H } _ { K }$ Inside each valley, the electrons are fully spin polarized due to the large spin-valley coupling, and opposite valleys have opposite spin due to time-reversal symmetry [30,31]. 

The parameters can be fixed by various approaches, among which the simplest one is to fit from DFT calculations for bilayer $\mathbf { M o T e } _ { 2 }$ at various stackings. This approach is adopted by Ref. [28] and the parameters are reproduced in Table I. The valence band structure with this set of parameters is shown in Fig. 1(b) at twist angle $3 . 8 9 ^ { \circ }$ The topmost valence band is dispersive with bandwidth larger than 20 meV. In addition, the two topmost bands overlap each other in energy. Both features are unfavorable for the emergence of FCIs. 

In this Letter, we seek to establish a better understanding of the single-particle band structure by performing largescale DFT calculations to take into account atomic relaxation, layer corrugation, and interlayer electric polarization (details in the Supplemental Material [32]). We choose the closest commensurate twist angle (3.89°) to the experimental value and construct the moir´e superlattice of MoTe using its monolayer unit cell with the optimized lattice constant $a = 3 . 5 2 \mathrm { ~ \AA ~ }$ The band structure of the moir´e superlattice is presented in Fig. 1(a) as red dots. The DFT result shows significant lattice reconstruction in both in plane and out of plane directions [Figs. 1(c) and 1(d)]. We then fit the continuum model parameters to the DFT band structure, and the result is presented in Table I. Compared to the parameters from Ref. [28], our parameters features a much larger moir´e potential and interlayer tunneling, which is likely caused by the significant lattice reconstruction [37–39], resulting in an isolated Chern band with bandwidth of roughly 9 meV [cf. Fig. 1(a)]. 


(a)



(b)


![](images/46d38716f7cf47bd38168263399bf3c3a7bde44c777ccf1f260269887383b97b.jpg)


![](images/71e3377192c725f36e6aabf509649acae964cb151fb169e1d9c19bab78e07ef6.jpg)



(c)



(d)


![](images/2cc00031f5557ef356942fa3b5ef91cfbeeda8fa269795fb8f58b985ed468def.jpg)


![](images/bc1aeaf55c77b7bc17df43082205beed4e918a0ff29cd77c527a6dbc05ce2b72.jpg)



FIG. 1. Band structures for K valley calculated by continuum model with parameters derived from our DFT calculation (a) and parameters from Ref. [28] (b). Chern numbers for the two topmost bands are labeled in the plot. Kohn-Sham DFT band structure is plotted with red circles in (a), and the two DFT bands labeled by empty circles are from the Γ valley. The twist angle is $3 . 8 9 ^ { \circ }$ . (c) In plane atomic displacement field in a moir´e unit cell after relaxation, and (d) shows interlayer (IL) distance for the line cut in (c). The maximal in plane atomic displacement is 0.32 Å. High symmetry stackings are labeled in (c). MM/XM/MX denotes the stacking where the metal/chalcogen/metal atoms of the top layer are directly above the metal/metal/chalcogen atoms of the bottom layer, respectively.


Fractional Chern insulator at $\nu = - 2 / 3$ .—Having established the existence of an isolated, relatively flat Chern band, we investigate whether FCIs can be stabilized. We adopt the following form of the Coulomb interaction: 

$$
H _ {\text { int }} = \frac {1}{2 A} \sum_ {l, l ^ {\prime}, \tau , \tau^ {\prime}, k, k ^ {\prime}, q} V (\boldsymbol {q}) c _ {l \tau k + q} ^ {\dagger} c _ {l ^ {\prime} \tau^ {\prime} k ^ {\prime} - q} ^ {\dagger} c _ {l ^ {\prime} \tau^ {\prime} k ^ {\prime}} c _ {l \tau k},\tag{2}
$$

where $V ( \pmb q ) = e ^ { 2 } \operatorname { t a n h } ( | q | d ) / 2 \epsilon _ { 0 } \epsilon | q |$ is the Coulomb interaction with dual-gate screening, A is the area of the system (proportional to the number of k points in the calculations), d is the distance between the twisted bilayer $\mathbf { M o T e } _ { 2 }$ and two symmetric metal gates, $\epsilon _ { 0 }$ is the vacuum permittivity, and ϵ is the relative dielectric constant. Here, $c _ { l \tau k } ^ { \dagger }$ creates a plane wave with momentum k at valley τ and layer l. 


(a)



(b)


![](images/9fce5d1ffea745025927f255404c1a12be99a27b706a7d5e9cbe770cd9a5d2bc.jpg)


![](images/f2e5dd9a3f9e06d1957f04e27012a896b74880c72d669e48f82ecdd0e72939a3.jpg)



(c)


![](images/8b70939f136dcbab17a1f620742886f7c6a6b61c3c4cc972ad6a8d2a9c735078.jpg)



(d)


![](images/f9574648e2be0c03b9b6e75e06fb75e2f1694d54f40baf9d81016e5205ae3b2f.jpg)



FIG. 2. For a large range of twist angles, the ground state of twisted bilayer $\mathbf { M o T e } _ { 2 }$ is fully valley polarized, as shown in the top of (a). The arrows represent spin, which is locked to the valley $( K$ and $K ^ { \prime } )$ degree of freedom. The bottom of (a) shows an example of the excited state with one spin flip. The energy difference between the lowest state with spin flip(s) and the ground state is defined as the spin gap, shown in (b) for $\nu = - 1 / 3$ and $\nu = - 2 / 3$ as a function of the twist angle. The calculation for (b) is carried out with $3 \times 4$ unit cells. (c) The many-body spectrum as a function of total crystalline momentum with the assumption of full valley polarization at $\nu = - 2 / 3$ . (d) The evolution of ground states under flux insertion along the $k _ { 2 }$ direction. During the flux insertion, the many-body gap is maintained. The calculation for (c) and (d) is carried out with $4 \times 6$ unit cells; the dielectric constant is chosen to be $^ { 1 5 ; }$ the distance between gate and sample is chosen to be $d = 3 0 0 \mathrm { ~ \AA } ;$ the twist angle is 3.89°.


Because of spin-valley locking, τ can also be understood as the spin label. We project the interaction onto the topmost moir´e band and carry out ED calculations. We choose $\epsilon = 1 5$ to make the characteristic interaction strength smaller than the averaged energy gap. For smaller ϵ, we present the phase diagram in the Supplemental Material [32]. While we left a more accurate treatment to include the band mixing for future studies, there are evidences that FCIs can still be stabilized even when the interaction exceeds the band gap [40,41]. 

The precursor to the FCI is spontaneous time-reversal symmetry breaking. We first perform ED calculations taking both valleys into account with a system size of $3 \times 4$ unit cells. We find that, over a broad range of twist angles, the ground state for both $\nu = - 1 / 3$ and $\nu = - 2 / 3$ is fully valley polarized, with holes occupying only one valley [cf. Fig. 2(a)]. Since the spin and valley indices are coupled, full valley polarization implies full spin polarization. The spin gap, defined as the energy difference between the lowest-energy state that does not exhibit full valley polarization and the fully valley-polarized ground state, is shown in Fig. 2(b) for both $\nu = - 1 / 3$ and $\nu = - 2 / 3$ . The spin gap for $\nu = - 1 / 3$ is much smaller than that of $\nu = - 2 / 3$ , indicating much weaker ferromagnetism of the former. The difference in the spin gap is consistent with the experimental observation that ferromagnetism appears at ∼4.5 K at $\nu = - 2 / 3$ , whereas no ferromagnetism is observed at $\nu = - 1 / 3$ down to base temperature of 1.6 K [42]. 

Given the large spin gap and strong ferromagnetism at $\nu = - 2 / 3$ , we further carry out ED calculations for a single valley, which allows us to consider a larger system with $4 \times 6$ unit cells. 

The most important signature of the FCI is the ground state degeneracy when the system is put on a torus [43,44]. Indeed, the ED energy spectrum shows three nearly degenerate states, separated by an energy gap from other states, as shown in Fig. 2(c). Under flux insertion, the three ground states evolve into each other, exhibiting a 6π periodicity [Fig. 2(d)]. We also calculate the many-body Chern number [43] at this filling to be $- 2 / 3$ , which is consistent with the experimental observation [27]. The single-particle occupation number, defined as $n ( k ) =$ $\langle c _ { k } ^ { \dagger } c _ { k } \rangle$ is presented in the Supplemental Material [32]. The uniformity of $n ( k )$ is a strong indicator favoring FCI over charge density wave states, one of the FCI’s competing phases. The above evidences provide strong evidence to the existence of the FCI in twisted $\mathbf { M o T e } _ { 2 }$ at $\nu = - 2 / 3$ 

Phase diagram.—The emergence of the FCI state depends on the dominance of electron-electron interaction energy over single-particle kinetic energy. This ratio between the two energy scales can be adjusted by two factors: the single-particle bandwidth and the environmental dielectric screening of electron-electron interactions. The most experimentally accessible knob to tune the bandwidth is changing the twist angle. For example, in twisted bilayer graphene, the flat band emerges at a series of magic angles [45]. For twisted bilayer TMD systems, the bandwidth is less sensitive to twist angles and we find relative isolated flat bands for twist angle $3 ^ { \circ } { \ - - } 4 ^ { \circ }$ . The screening of the electron-electron interactions can be tuned by changing ϵ, as well as the distance between the sample and the metal gate d: larger d leads to weaker screening of the electron-electron interactions. 

The phase diagram, as a function of θ and either ϵ or d, is presented in Fig. 3. Within the range of the twist angle presented in Fig. 3, the bandwidth of the topmost valence band increases monotonically with θ. There are three distinct phases: the FCI phase, the valley-polarized (VP) phase, and the non-valley-polarized (NVP) phase. The NVP phase emerges in the weak electron-electron interaction regime, characterized by large θ and large ϵ [see Fig. 3(a)], or small d [see Fig. 3(b)]. For stronger interactions, all holes occupy the same valley, which is shown as the VP phase in Fig. 3. Our numerical evidence suggests that the VP states are most likely metal states with fully polarized spin, meaning they represent a half-metal state [23]. The nature of these VP phases is left for further studies. The FCI emerges from the VP phases at even stronger interactions. The many-body gap for the FCI phase has a peak at $\theta = 3 . 5 ^ { \circ }$ , indicating an optimal twist angle for the observation of FCIs. This optimal angle is close to the twist angle (∼3.7°) of the device in which the FCI is observed [27]. 

![](images/9279fd47c35e61ada960191556fd54bf45cb3020e5f38c2e5529b45de5d7607a.jpg)


![](images/b8abfa81b2f9546906d3faf6fb1ae34db4f51dc1a398e9ea0f51f78c125e3ca7.jpg)



FIG. 3. (a) Phase diagram as a function of twist angle and dielectric constant at $\nu = - 2 / 3$ . d is chosen as 300 Å. (b) Phase diagram as a function of twist angle and d at $\nu = - 2 / 3$ . ϵ is chosen as 15. Many-body gap is shown in both (a) and (b) by color in the FCI phase (units in meV). The FCI is identified with $4 \times 6$ unit cells. Valley polarization is identified with $3 \times 4$ unit cells.


The effect ofelectricfield.—An out of plane electric field generates potential differences between the top and bottom layers. Experimentally, it is observed that the FCI at $\nu = - 2 / 3$ can be suppressed by this out of plane electric field. This observation is not necessarily surprising, since layer potential differences will induce a topological phase transition at the single-particle level, making the topmost valence band topologically trivial. However, our calculations show that the FCI is suppressed well before the single-particle topological phase transition. In Fig. 4(a), we show that the many-body gap closes at electric field $E = 1 . 2 6 ~ \mathrm { m e V / \mathring { A } }$ In Fig. 4(b), we present the singleparticle band structure for $E = 1 . 3 9 ~ \mathrm { m e V / \mathring { A } }$ , where an isolated flat Chern band can still be observed. The bandwidth of the Chern band at $E = 1 . 3 9 ~ \mathrm { m e V / \mathring { A } }$ is comparable to that at $E = 0 . 0 ~ \mathrm { m e V / \mathring { A } }$ , but the FCI is already destroyed. In experiment, the ferromagnetism disappears at $E \sim 5 \mathrm { \ m e V / \mathring { A } } [ 2 7 ]$ , implying the FCI is destroyed at a smaller electric field. Our critical electric field $E =$ $1 . 2 6 ~ \mathrm { m e V / \mathring { A } }$ is consistent with this observation. 

It is well known that completely quenched kinetic energy (i.e., vanishing bandwidth) does not ensure the existence of the FCI, and a number of proposals [46–50] are put forward to identify the conditions for the FCI to emerge. Many of these proposals aim to design wave functions in a flat Chern band such that they closely resemble the wave functions of a Landau level. For example, since the Berry curvature Ω and quantum metric tensor g are constant for Landau levels, the flatness of these two quantities [46–48] in the reciprocal space is heuristically viewed as a promising indicator for the emergence of FCIs. In Fig. 4(c), we present the distributions of Ω and $\operatorname { t r } ( g )$ for $E = 0 . 0$ and $E = 1 . 3 9 ~ \mathrm { m e V / \mathring { A } }$ A serious deterioration of the flatness of the Berry curvature and quantum metric tensor can be observed when the FCI is destroyed. Both Ω and $\operatorname { t r } ( g )$ at $E = 1 . 3 9 ~ \mathrm { m e V / \mathring { A } }$ are concentrated at the moir´e K valley [Fig. 4(c)], which explains the suppression of the FCI state. 

![](images/abe302728c19f012ac986a026b74c04062fd6d79153abc98f1c998f5c2238cd3.jpg)


![](images/506c99bf7ae45b3b1c87ef5b1c69da2308607df26f967793691b707a902dca4a.jpg)


![](images/4ce57e993812edb56c2c5b3ef2b2d78499ca53126a9c25077921e993ec2b3f64.jpg)



FIG. 4. (a) Many-body gap as a function of external electric field at $\epsilon = 1 5$ and $\theta = 3 . 8 9 ^ { \circ }$ . (b) Band structure at $E =$ $1 . 3 9 \mathrm { \ m V / \mathring A }$ and $\theta = 3 . 8 9 ^ { \circ }$ (the FCI is already destroyed by the electric field). (c) Indicators for the emergence of the FCI at $E = 0 ~ \mathrm { m V / \mathring { A } }$ (left) and $E = 1 . 3 9 ~ \mathrm { m V / \mathring { A } }$ (right). Ω is the Berry curvature and trðgÞ is the trace of quantum metric tensor. The unit for (c) is the inverse of the area of the moir´e Brillouin zone.


Fractional Chern insulator at $\nu = - 3 / 5$ .—In addition to the $\nu = - 1 / 3$ state, the experimental observation also includes the $\nu = - 3 / 5$ state [27]. However, in our calculations, we do not observe a clear many-body gap at $\nu = - 3 / 5$ for the dielectric constant of $\epsilon = 1 5$ Nonetheless, when the dielectric constant is increased to $\epsilon = 8 .$ , which is closer to the experimental value, we find evidence for the FCI [32]. The characteristic interaction at $\epsilon = 8$ is larger than the energy gap between the two topmost bands, but the result can still be viewed as supporting evidences of the existence of the FCI at $\nu = - 3 / 5$ We also performed ED calculations for $\nu = - 1 / 5 , \ - 2 / 5$ , and $- 4 / 5$ and find no clear evidence of FCIs for $\nu = - 1 / 5$ and $- 4 / 5$ . For $\nu = - 2 / 5$ , we find evidence of a FCI phase with the many-body energy gap roughly half of that at $\nu = - 3 / 5$ 

In summary, we have presented a comprehensive theoretical study of the recent experimental observation of a fractional Chern insulator at zero magnetic field in a twisted $\mathbf { M o T e } _ { 2 }$ bilayer. Density-functional theory calculations reveal the existence of an isolated flat Chern band, which allows us to confirm the existence of a $\nu = - 2 / 3$ fractional Chern insulator state using exact diagonalization. Phase diagrams are presented to guide future experimental study of fractional Chern insulators in this system, and the suppression of fractional Chern insulators by external electric field is studied. Our findings offer valuable insights into the nature and properties of fractional Chern insulators in moir´e superlattices. 

Note added.—We recently became aware of an independent work on similar topics [51]. A recent ED calculation [52] based on our parameters has highlighted the importance of remote band effects. Another recent DFT calculation [53] using a different package has produced similar band structures as ours. 

We thank Jie Wang, Jiaqi Cai, and Heqiu Li for stimulating discussions. The exact diagonalization study is supported by DOE Award No. DE-SC0012509. The density-functional theory calculation is supported by the Center on Programmable Quantum Materials, an Energy Frontier Research Center funded by DOE BES under Award No. DE-SC0019443. Y. H. acknowledges support from the European Research Council (ERC) under the European Union Horizon 2020 Research and Innovation Programme (Grant Agreement No. 804213-TMCS). This work was facilitated through the use of advanced computational, storage, and networking infrastructure provided by the Hyak supercomputer system and funded by the University of Washington Molecular Engineering Materials Center at the University of Washington (DMR-1719797). 



tingcao@uw.edu <sup>†</sup>dixiao@uw.edu 





[1] D. C. Tsui, H. L. Stormer, and A. C. Gossard, Twodimensional magnetotransport in the extreme quantum limit, Phys. Rev. Lett. 48, 1559 (1982). 





[2] R. B. Laughlin, Anomalous quantum Hall effect: An incompressible quantum fluid with fractionally charged excitations, Phys. Rev. Lett. 50, 1395 (1983). 





[3] S. A. Parameswaran, R. Roy, and S. L. Sondhi, Fractional quantum Hall physics in topological flat bands, C.R. Phys. 14, 816 (2013). 





[4] E. J. Bergholtz and Z. Liu, Topological flat band models and fractional Chern insulators, Int. J. Mod. Phys. B 27, 1330017 (2013). 





[5] T. Neupert, C. Chamon, T. Iadecola, L. H. Santos, and C. Mudry, Fractional (Chern and topological) insulators, Phys. Scr. 2015, 014005 (2015). 





[6] Z. Liu and E. J. Bergholtz, Recent developments in fractional Chern insulators, in Reference Module in Materials Science and Materials Engineering (Elsevier, New York, 2023), p. B9780323908009001360. 





[7] S. Das Sarma, M. Freedman, and C. Nayak, Topologically protected qubits from a possible non-Abelian fractional quantum Hall state, Phys. Rev. Lett. 94, 166802 (2005). 





[8] C. Nayak, S. H. Simon, A. Stern, M. Freedman, and S. Das Sarma, Non-Abelian anyons and topological quantum computation, Rev. Mod. Phys. 80, 1083 (2008). 





[9] T. Neupert, L. Santos, C. Chamon, and C. Mudry, Fractional quantum Hall states at zero magnetic field, Phys. Rev. Lett. 106, 236804 (2011). 





[10] E. Tang, J.-W. Mei, and X.-G. Wen, High-temperature fractional quantum Hall states, Phys. Rev. Lett. 106, 236802 (2011). 





[11] D. Sheng, Z.-C. Gu, K. Sun, and L. Sheng, Fractional quantum Hall effect in the absence of Landau levels, Nat. Commun. 2, 389 (2011). 





[12] N. Regnault and B. A. Bernevig, Fractional Chern insulator, Phys. Rev. X 1, 021014 (2011). 





[13] D. Xiao, W. Zhu, Y. Ran, N. Nagaosa, and S. Okamoto, Interface engineering of quantum Hall effects in digital transition metal oxide heterostructures, Nat. Commun. 2, 596 (2011). 





[14] E. Y. Andrei, D. K. Efetov, P. Jarillo-Herrero, A. H. MacDonald, K. F. Mak, T. Senthil, E. Tutuc, A. Yazdani, and A. F. Young, The marvels of moir´e materials, Nat. Rev. Mater. 6, 201 (2021). 





[15] D. M. Kennes, M. Claassen, L. Xian, A. Georges, A. J. Millis, J. Hone, C. R. Dean, D. Basov, A. N. Pasupathy, and A. Rubio, Moir´e heterostructures as a condensed-matter quantum simulator, Nat. Phys. 17, 155 (2021). 





[16] P. J. Ledwith, G. Tarnopolsky, E. Khalaf, and A. Vishwanath, Fractional Chern insulator states in twisted bilayer graphene: An analytical approach, Phys. Rev. Res. 2, 023237 (2020). 





[17] C. Repellin and T. Senthil, Chern bands of twisted bilayer graphene: Fractional Chern insulators and spin phase transition, Phys. Rev. Res. 2, 023238 (2020). 





[18] A. Abouelkomsan, Z. Liu, and E. J. Bergholtz, Particle-hole duality, emergent fermi liquids, and fractional Chern insulators in moir´e flatbands, Phys. Rev. Lett. 124, 106803 (2020). 





[19] P. Wilhelm, T. C. Lang, and A. M. Läuchli, Interplay of fractional Chern insulator and charge density wave phases in twisted bilayer graphene, Phys. Rev. B 103, 125406 (2021). 





[20] Z. Liu, A. Abouelkomsan, and E. J. Bergholtz, Gate-tunable fractional Chern insulators in twisted double bilayer graphene, Phys. Rev. Lett. 126, 026801 (2021). 





[21] Y. Sheffer and A. Stern, Chiral magic-angle twisted bilayer graphene in a magnetic field: Landau level correspondence, exact wave functions, and fractional Chern insulators, Phys. Rev. B 104, L121405 (2021). 





[22] H. Li, U. Kumar, K. Sun, and S.-Z. Lin, Spontaneous fractional Chern insulators in transition metal dichalcogenide moir´e superlattices, Phys. Rev. Res. 3, L032070 (2021). 





[23] V. Cr´epel and L. Fu, Anomalous Hall metal and fractional Chern insulator in twisted transition metal dichalcogenides, Phys. Rev. B 107, L201109 (2023). 





[24] N. Morales-Durán, J. Wang, G. R. Schleder, M. Angeli, Z. Zhu, E. Kaxiras, C. Repellin, and J. Cano, Pressureenhanced fractional Chern insulators in moir´e transition metal dichalcogenides along a magic line, Phys. Rev. Res. 5, 032022 (2023). 





[25] Y. Xie, A. T. Pierce, J. M. Park, D. E. Parker, E. Khalaf, P. Ledwith, Y. Cao, S. H. Lee, S. Chen, P. R. Forrester et al., Fractional Chern insulators in magic-angle twisted bilayer graphene, Nature (London) 600, 439 (2021). 





[26] E. M. Spanton, A. A. Zibrov, H. Zhou, T. Taniguchi, K. Watanabe, M. P. Zaletel, and A. F. Young, Observation of fractional Chern insulators in a Van der Waals heterostructure, Science 360, 62 (2018). 





[27] J. Cai, E. Anderson, C. Wang, X. Zhang, X. Liu, W. Holtzmann, Y. Zhang, F. Fan, T. Taniguchi, K. Watanabe, Y. Ran, T. Cao, L. Fu, D. Xiao, W. Yao, and X. Xu, Signatures of fractional quantum anomalous Hall states in twisted $\begin{array} { r } { \mathbf { M o T e } _ { 2 } , } \end{array}$ , Nature (London) 622, 63 (2023). 





[28] F. Wu, T. Lovorn, E. Tutuc, I. Martin, and A. H. MacDonald, Topological insulators in twisted transition metal dichalcogenide homobilayers, Phys. Rev. Lett. 122, 086402 (2019). 





[29] H. Yu, M. Chen, and W. Yao, Giant magnetic field from moir´e induced berry phase in homobilayer semiconductors, Natl. Sci. Rev. 7, 12 (2020). 





[30] D. Xiao, G.-B. Liu, W. Feng, X. Xu, and W. Yao, Coupled spin and valley physics in monolayers of mos and other group-VI dichalcogenides, Phys. Rev. Lett. 108, 196802 (2012). 





[31] T. Cao, G. Wang, W. Han, H. Ye, C. Zhu, J. Shi, Q. Niu, P. Tan, E. Wang, B. Liu et al., Valley-selective circular dichroism of monolayer molybdenum disulphide, Nat. Commun. 3, 887 (2012). 





[32] See Supplemental Material at http://link.aps.org/ supplemental/10.1103/PhysRevLett.132.036501 for DFT details and additional calculations, which includes Refs. [1,33–35]. 





[33] J. M. Soler, E. Artacho, J. D. Gale, A. García, J. Junquera, P. Ordejón, and D. Sánchez-Portal, The SIESTA method for ab initio order-N materials simulation, J. Condens. Matter Phys. 14, 2745 (2002). 





[34] D. R. Hamann, Optimized norm-conserving Vanderbilt pseudopotentials, Phys. Rev. B 88, 085117 (2013). 





[35] J. P. Perdew, K. Burke, and M. Ernzerhof, Generalized gradient approximation made simple, Phys. Rev. Lett. 77, 3865 (1996). 





[36] S. Grimme, Semiempirical GGA-type density functional constructed with a long-range dispersion correction, J. Comput. Chem. 27, 1787 (2006). 





[37] M. H. Naik and M. Jain, Ultraflatbands and shear solitons in moir´e patterns of twisted bilayer transition metal dichalcogenides, Phys. Rev. Lett. 121, 266401 (2018). 





[38] E. Li, J.-X. Hu, X. Feng, Z. Zhou, L. An, K. T. Law, N. Wang, and N. Lin, Lattice reconstruction induced multiple ultra-flat bands in twisted bilayer $\mathrm { W S e } _ { 2 } ,$ Nat. Commun. 12, 5601 (2021). 





[39] H. Li, S. Li, M. H. Naik, J. Xie, X. Li, J. Wang, E. Regan, D. Wang, W. Zhao, S. Zhao et al., Imaging moir´e flat bands in three-dimensional reconstructed WSe =WS superlattices, Nat. Mater. 20, 945 (2021). 





[40] S. Kourtis, T. Neupert, C. Chamon, and C. Mudry, Fractional Chern insulators with strong interactions that far exceed band gaps, Phys. Rev. Lett. 112, 126806 (2014). 





[41] A. G. Grushin, J. Motruk, M. P. Zaletel, and F. Pollmann, Characterization and stability of a fermionic $\nu = 1 / 3$ fractional Chern insulator, Phys. Rev. B 91, 035136 (2015). 





[42] Jiaqi Cai and Xiaodong Xu (private communications). 





[43] Q. Niu, D. J. Thouless, and Y.-S. Wu, Quantized Hall conductance as a topological invariant, Phys. Rev. B 31, 3372 (1985). 





[44] X. G. Wen and Q. Niu, Ground-state degeneracy of the fractional quantum Hall states in the presence of a random potential and on high-genus Riemann surfaces, Phys. Rev. B 41, 9377 (1990). 





[45] R. Bistritzer and A. H. MacDonald, Moir´e bands in twisted double-layer graphene, Proc. Natl. Acad. Sci. U.S.A. 108, 12233 (2011). 





[46] S. H. Simon and M. S. Rudner, Contrasting lattice geometry dependent versus independent quantities: Ramifications for berry curvature, energy gaps, and dynamics, Phys. Rev. B 102, 165148 (2020). 





[47] S. A. Parameswaran, R. Roy, and S. L. Sondhi, Fractional Chern insulators and the $W _ { \infty }$ algebra, Phys. Rev. B 85, 241308(R) (2012). 





[48] R. Roy, Band geometry of fractional topological insulators, Phys. Rev. B 90, 165139 (2014). 





[49] J. Wang, J. Cano, A. J. Millis, Z. Liu, and B. Yang, Exact Landau level description of geometry and interaction in a flatband, Phys. Rev. Lett. 127, 246403 (2021). 





[50] P. J. Ledwith, A. Vishwanath, and D. E. Parker, Vortexability: A unifying criterion for ideal fractional Chern insulators, Phys. Rev. B 108, 205144 (2023). 





[51] A. Reddy, F. Alsallom, Y. Zhang, T. Devakul, and L. Fu, Fractional quantum anomalous Hall states in twisted bilayer MoTe<sub>2</sub> and WSe<sub>2</sub>, Phys. Rev. B 108, 085117 (2023). 





[52] J. Yu, J. Herzog-Arbeitman, M. Wang, O. Vafek, B. A. Bernevig, and N. Regnault, Fractional Chern insulators vs non-magnetic states in twisted bilayer MoTe<sub>2</sub>, arXiv: 2309.14429. 





[53] Y. Jia, J. Yu, J. Liu, J. Herzog-Arbeitman, Z. Qi, N. Regnault, H. Weng, B. A. Bernevig, and Q. Wu, Moir´e fractional Chern insulators I: First-principles calculations and continuum models of twisted bilayer MoTe<sub>2</sub>, arXiv: 2311.04958. 

