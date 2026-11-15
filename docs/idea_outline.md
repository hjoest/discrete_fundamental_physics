# Idea Outline: Discrete Fundamental Physics from Posets and Finite Groups

## 1. Foundational Premise (Strict Finitism & Information Simplicity)

A fundamental physical theory should rest on the smallest, simplest set of axioms.
- **Simple primitives:** Countable sets, preorders / partially ordered sets (posets), and finite discrete groups.
- **Non-simple constructs (emergent only):** Continua ($\mathbb{R}^n$), smooth manifolds, continuous Lie groups, and infinite-dimensional vector spaces. A single real coordinate carries infinite Shannon information; therefore, continuum structures must arise only as large-scale statistical or coarse-grained approximations, never as microscopic axioms.

---

## 2. Kinematic Structure: Causal Base + Discrete Internal Fiber

The kinematic arena combines two discrete structures:
1. **Causal Poset $(P, \le_P)$:** A countable, locally finite partially ordered set whose elements represent spacetime events and whose partial order encodes irreversible causal precedence (as in Sorkin's Causal Set Theory, *"Order + Number = Geometry"*).
2. **Finite Group $G$:** A finite non-abelian group encoding internal degrees of freedom (spin, chirality, flavor, and gauge charges).

### Event Space $Z$ and Fiber Incomparability
To unify spacetime causality with internal gauge structure in the spirit of Kaluza–Klein theory and Connes' Noncommutative Geometry ($M \times F \to P \times G$) without introducing instantaneous closed timelike curves (CTCs):
- Events live in $Z \subseteq P \times G$ with natural projections $\pi_P: Z \to P$ and $\pi_G: Z \to G$.
- **Strict Causal Inheritance (Incomparable Fibers):** Distinct states in the same fiber $F_p = \pi_P^{-1}(p)$ are mutually incomparable (spacelike):
  $$z_1 <_Z z_2 \implies \pi_P(z_1) <_P \pi_P(z_2).$$
  Together with reflexivity $z \le_Z z$, this guarantees that $(Z, \le_Z)$ is a strict **partially ordered set** (antisymmetric, transitive, locally finite).

### Physical Dictionary of the Product Structure
| Geometric Component in $(Z, \prec_Z)$ | Physical Counterpart in Emergent QFT / Standard Model |
| --- | --- |
| **Base poset $(P, \prec_P)$** | 4D Lorentzian spacetime causal structure & metric volume |
| **Covering links along $P$** | Spacetime propagation & Yang–Mills gauge connections $A_\mu$ |
| **Fiber coordinates in $F_p \subseteq \{p\} \times G$** | Fermion & boson multiplets evaluated at internal sites |
| **Internal transitions across $G$** | Discrete Higgs field $\Phi$ & Yukawa mass couplings (Connes–Lott mechanism) |

---

## 3. Microscopic Dynamics: The Ordered Product Conservation Law

For each event $z \in Z$, let its immediate predecessors (incoming covering links) and immediate successors (outgoing covering links) be:
$$\mathrm{Pred}(z) = \{u \in Z \mid u \prec_Z z\} = \{u_1, \dots, u_m\}, \qquad \mathrm{Succ}(z) = \{v \in Z \mid z \prec_Z v\} = \{v_1, \dots, v_n\}.$$

**Axiom (Local Algebraic Balance):** At every interior event $z \in Z$, there exist orderings (permutations $\sigma \in S_m$ and $\tau \in S_n$) such that the ordered product of incoming group labels equals the ordered product of outgoing group labels in $G$:
$$\prod_{i=1}^m \pi_G(u_{\sigma(i)}) = \prod_{j=1}^n \pi_G(v_{\tau(j)}).$$

### Consequences of the Conservation Law
1. **Abelian Sector (Exact Noether / Kirchhoff Law):** Projecting to the abelianization $G^{\text{ab}} = G/[G,G]$ eliminates ordering dependence, yielding an exact discrete continuity equation $\sum_{u \prec z} [\pi_G(u)] = \sum_{z \prec v} [\pi_G(v)]$ for conserved abelian charges.
2. **Non-Abelian Sector (Discrete Yang–Mills / Intertwiner Constraint):** When $G$ is non-abelian, the sensitivity to factor ordering is governed by commutators $[g_i, g_j] \in [G,G]$, generating non-linear self-interactions analogous to non-abelian Gauss laws and spin-foam vertex intertwiners.
3. **Permutation Gauge Freedom:** The existence of *some* ordering $(\sigma, \tau) \in S_m \times S_n$ rather than a rigid planar embedding treats incident-edge permutations as a local discrete gauge symmetry at each vertex.
4. **Dynamic Poset Filter (Suppressing Kleitman–Rothschild Pathologies):** Requiring valid non-trivial group valuations acts as a combinatorial constraint on admissible posets $P$, filtering out pathological causal sets and favoring regular, high-connectivity aperiodic/quasicrystalline graphs.

---

## 4. Candidate Finite Groups $G$ and Falsifiability

To recover Alain Connes' finite spectral triple algebra for the Standard Model,
$$\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C}),$$
the finite group $G$ must generate the required real/complex/quaternionic matrix algebras (via its real group algebra $\mathbb{R}[G]$) and fermion representations:

| Candidate Group $G$ | Order $\lvert G \rvert$ | Role & Motivation |
| --- | --- | --- |
| **$2I \cong \mathrm{SL}(2, \mathbb{F}_5)$** | $120$ | Binary icosahedral group (vertices of the regular 4D 600-cell); maximal finite subgroup of $\mathrm{SU}(2) \cong U(1, \mathbb{H})$; generates quaternionic ($\mathbb{H}$) electroweak doublets. |
| **$\mathrm{PSL}(2, \mathbb{F}_7) \cong \mathrm{GL}(3, \mathbb{F}_2)$** | $168$ | Klein quartic / Fano plane automorphism group; discrete symmetry of octonion multiplication ($\mathbb{O}$), linked to $\mathrm{SU}(3)_c$ color. |
| **$\Delta(27) \cong (\mathbb{Z}_3 \times \mathbb{Z}_3) \rtimes \mathbb{Z}_3$** | $27$ | Discrete Heisenberg group ($SU(3)$ subgroup); smallest group with faithful 3D complex irreps and triality. |
| **$\Sigma(216)$ (Hessian group)** | $216$ (or $648$ in $\mathrm{SU}(3)$) | Principal finite subgroup of $\mathrm{PSU}(3)$ / $\mathrm{SU}(3)$ used in discrete flavor/color models. |
| **$2I \times \mathrm{PSL}(2, \mathbb{F}_7)$** | $20,160$ | Minimal bit-depth composite candidate combining quaternionic weak isospin ($2I$) and octonionic color ($\mathrm{PSL}(2, \mathbb{F}_7)$). |
| **$2I \times \Sigma(216)$ / $2I \times \Delta(27)$** | $25,920$ / $3,240$ | Direct geometric product of binary icosahedral weak sector with discrete $\mathrm{SU}(3)$ color/generation sector. |
| **$W(E_8)$** | $696,729,600$ | Weyl group of $E_8$; single unified reflection group preserving the $E_8$ root lattice, containing $W(\mathrm{Spin}(10)) \times \mathbb{Z}_3$ and the H4 (600-cell) symmetries. |

### Quantitative Falsification Targets
Because finite groups possess no continuous free parameters, candidate groups $G$ can be tested and ruled out by computing:
1. **Real Wedderburn Decomposition of $\mathbb{R}[G]$:** Checking whether $\mathbb{R}[G]$ contains $\mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ (using Frobenius–Schur indicators $\iota \in \{+1, 0, -1\}$) and counting the multiplicity of fundamental irreps (number of fermion generations $N_g = 3$).
2. **Weinberg Angle $\sin^2\theta_W(\Lambda)$ & Coupling Ratios:** Evaluating generator trace normalizations $\mathrm{Tr}(T_{3L}^2)/\mathrm{Tr}(Q^2)$ across the fiber Hilbert space.
3. **Bare Higgs Quartic Ratio from Cayley Spectra:** Computing the eigenvalue spectrum of the discrete Dirac / Laplacian operator $D_G$ on the Cayley graph of $(G, S)$ to evaluate $\kappa_G = \mathrm{Tr}(D_G^4) / (\mathrm{Tr}(D_G^2))^2$ and the resulting bare Higgs mass prediction.
4. **Combinatorial Branching & Mass Hierarchies:** Analyzing the solution space of $\prod \pi_G(u_{\sigma(i)}) = \prod \pi_G(v_{\tau(j)})$ across conjugacy classes and word-length shells of $G$.
