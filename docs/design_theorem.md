# The 4-Point Design Theorem for Viable Finite Groups in $(P, G, Z, \prec_Z)$

This document formalizes the **4-Point Design Theorem** derived from the seven falsification tests in [`falsification_report.md`](falsification_report.md) that eliminated both Option A ($\Sigma(216 \times 3)$) and Option B ($2I \times \Sigma(216 \times 3)$).

---

## 1. Summary of the Four Necessary Conditions

For a discrete spacetime-matter tuple $(P, G, Z, \prec_Z)$ governed by the ordered product conservation law and the real group algebra $\mathbb{R}[G]$ to simultaneously reproduce:
1. **4D aperiodic spacetime geometry and chiral Lorentz Weyl spinors**,
2. **Connes' Standard Model gauge algebra $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ with anomaly cancellation and $\sin^2\theta_W(\Lambda) = 3/8$**,
3. **Three universal generations ($N_g = 3$) of both quarks and leptons over $\mathbb{R}$ with non-trivial CKM/PMNS mixing ($V_{\text{CKM}} \neq I_3$)**, and
4. **Non-trivial Kirchhoff charge conservation and resistance to Kleitman–Rothschild causal washout**,

the finite group $G$ must satisfy four exact structural conditions:

| Condition | Mathematical Requirement on $G$ | What Fails If Violated (No-Go Theorem) |
| :--- | :--- | :--- |
| **1. External 4D Aperiodic Spinor Factor** | $G$ must contain a dedicated spacetime factor $G_{\text{st}} \cong 2I \subset \mathrm{Spin}(4)$ with two faithful quaternionic doublets $\mathbf{2}_{\text{spin}, L} \oplus \mathbf{2}_{\text{spin}, R}$ ($\iota = -1$) having irrational traces in $\mathbb{Q}(\sqrt{5})$ whose direct sum embeds into the 8D $E_8$ lattice ($\mathrm{Tr}_{\mathbb{R}^8} \in \mathbb{Z}$), and **neither doublet may be reused as weak isospin**. | **4D Crystallographic Trap & Spin–Isospin Dilemma (Tests 1 & 7):** Without irrational 4D traces, the cut-and-project tiling collapses to a periodic crystal (e.g., $2T \to F_4$ 24-cell lattice). Reusing $\mathbf{2}_{\text{spin}, L}$ for $\mathrm{SU}(2)_{W,L}$ locks weak isospin to spacetime chirality. |
| **2. Non-Abelian 3D Flavor Triplication** | The 3 fermion generations ($N_g = 3$) must transform as a **faithful irreducible 3-dimensional complex representation $\mathbf{3}_f$ ($\iota(\mathbf{3}_f) = 0$)** of a non-abelian flavor group $G_{\text{flavor}}$ with a non-trivial $\mathbb{Z}_3$ center $Z(G_{\text{flavor}}) \cong \mathbb{Z}_3$ (minimal order: $\lvert G_{\text{flavor}}\rvert = 27$ by Schur's center-index theorem, realized with exponent 3 by $\Delta(27)_f$), **never** as an orbit $\rho^{(0)} \otimes \lambda^a$ of a 1D abelian character $\lambda \in G^{\text{ab}}$. | **Lepton Collapse & Circulant CKM No-Go Theorems (Tests 5 & 6):** Any 1D $\mathbb{Z}_3$-orbit forces $\overline{\rho_{\text{singlet}} \otimes \lambda} = \rho_{\text{singlet}} \otimes \lambda^2$ (collapsing 3 lepton generations into $1$ real + $1$ complex over $\mathbb{R}$) and forces $M_u, M_d$ to be $\mathbb{Z}_3$-circulant ($\lvert V_{\text{CKM}}\rvert = I_3$). |
| **3. Complete Internal Chiral Gauge Factor** | The internal gauge factor $G_{\text{gauge}}$ must independently provide $\mathrm{SU}(2)_{W,L} \times \mathrm{SU}(2)_{W,R} \times \mathrm{SU}(3)_c$ (or distinct $u_R, d_R$ triplets alongside $\mathrm{SU}(2)_{W,L}$) so that $N(\mathbf{3}_L) = N(\mathbf{3}_R) = 6$, all chiral gauge anomalies cancel internally, and $\sin^2\theta_W(\Lambda) = 3/8$. | **Missing Right-Handed Quark & $\mathrm{SU}(3)_c^3$ Anomaly Trap (Tests 2 & 3):** A single internal $\mathrm{SU}(2)$ with one $\mathbf{3}_c$ per generation leaves half of $u_R/d_R$ missing and yields $\mathcal{A}(\mathrm{SU}(3)_c^3) = 6 - 3 = +3 \neq 0$. |
| **4. Non-Perfect Abelianization & Shared-Permutation Locking** | $G$ must have a non-trivial abelianization $G^{\text{ab}} = G/[G,G] \neq \{1\}$ and a commutator coset size $\lvert [G,G] \rvert \gg d!$ for $d \in \{3, 4\}$. | **Kleitman–Rothschild Causal Washout (Test 4):** Perfect groups ($G^{\text{ab}} = \{1\}$) have zero conserved Kirchhoff charge and suffer $100\%$ causal washout ($S_d \to 0\text{ bits}$) at degree $d \ge 5$. |

---

## 2. Mathematical Proofs of the Four Design Conditions

### Theorem 2.1 (Condition 1: Irrational 4D Cut-and-Project & Spin–Isospin Independence)
Let $\rho_{\parallel}: G \to \mathrm{O}(4)$ be the physical spacetime projection of a cut-and-project scheme from a periodic lattice $\Lambda \subset \mathbb{R}^D \cong \mathbb{R}^4_{\parallel} \oplus \mathbb{R}^{D-4}_{\perp}$.
1. **Aperiodicity Requirement:** By the 4D Crystallographic Restriction Theorem, if $\mathrm{Tr}(\rho_{\parallel}(g)) \in \mathbb{Z}$ for all $g \in G$, then $\rho_{\parallel}(G)$ preserves a rank-4 periodic lattice in $\mathbb{R}^4_{\parallel}$. To obtain a strictly aperiodic 4D quasicrystal with minimal embedding dimension $D = 8$, $\rho_{\parallel}$ must have quadratic irrational traces in $\mathbb{Q}(\sqrt{5})$ paired with its Galois conjugate $\rho_{\perp}$ so that $\mathrm{Tr}(\rho_{\parallel}(g)) + \mathrm{Tr}(\rho_{\perp}(g)) \in \mathbb{Z}$. The unique maximal finite subgroup of $\mathrm{Spin}(4) \cong \mathrm{SU}(2)_{\text{spin}, L} \times \mathrm{SU}(2)_{\text{spin}, R}$ realizing this in each chiral factor is the binary icosahedral group $2I$ (order 120), via $\mathbf{2}_{\text{spin}, L} \oplus \mathbf{2}_{\text{spin}, R} \cong \mathbb{R}^4_{\parallel} \oplus \mathbb{R}^4_{\perp} \supset \Lambda_{E_8}$.
2. **Spin–Isospin Independence:** Because left-handed Standard Model fermions include both weak doublets ($Q_L, L_L$) and weak singlets ($u_R^c, d_R^c, e_R^c, \nu_R^c$ in left-handed Weyl notation), the Lorentz spinor group $\mathrm{Spin}(4)$ must commute with internal weak isospin $\mathrm{SU}(2)_{W,L}$. Therefore, $G_{\text{st}} \cong 2I$ cannot be reused as $\mathrm{SU}(2)_{W,L}$. $\blacksquare$

---

### Theorem 2.2 (Condition 2: Why Generations Must Form an Irreducible Non-Abelian Triplet $\mathbf{3}_f$)
Suppose a candidate group attempts to generate $N_g = 3$ generations by tensoring a single-generation representation $\rho^{(0)}$ with the 1-dimensional characters $\{1, \lambda, \lambda^2\}$ of a $\mathbb{Z}_3$ quotient ($\rho^{(a)} = \rho^{(0)} \otimes \lambda^a$, as in $\Sigma(216 \times 3)$). Then:
1. **Lepton Collapse over $\mathbb{R}$:** For any color-singlet representation $\rho_0$ that is self-conjugate ($\bar{\rho}_0 \cong \rho_0$, $\iota(\rho_0) = \pm 1$), complex conjugation swaps $a=1$ and $a=2$:
   $$\overline{\rho_0 \otimes \lambda} = \bar{\rho}_0 \otimes \bar{\lambda} \cong \rho_0 \otimes \lambda^2.$$
   In $\mathbb{R}[G]$, $(\rho_0 \otimes \lambda) \oplus (\rho_0 \otimes \lambda^2)$ forms a **single** complex matrix block $M_d(\mathbb{C})$ (one Dirac generation + its antiparticle), while $\rho_0 \otimes 1$ remains real/quaternionic ($M_d(\mathbb{R})$ or $M_{d/2}(\mathbb{H})$). Thus leptons fail to form 3 complex generations.
2. **Circulant CKM Triviality:** Because $\overline{\rho_L^{(a)}} \otimes \rho_R^{(b)} = (\overline{\rho_L^{(0)}} \otimes \rho_R^{(0)}) \otimes \lambda^{b-a}$, every $G$-invariant Yukawa matrix $M_{ab}$ depends only on $(b - a) \bmod 3$. Since all $3 \times 3$ circulant matrices commute and share the Fourier eigenvector matrix $U_\omega$, $V_{\text{CKM}} = U_u^\dagger U_d = I_3$.

**Resolution via Non-Abelian $\mathbf{3}_f$ ($\iota(\mathbf{3}_f) = 0$):**
Let the 3 generations transform as a faithful 3D complex irrep $\mathbf{3}_f$ ($\iota(\mathbf{3}_f) = 0$) of the discrete Heisenberg group $G_{\text{flavor}} = \Delta(27)_f \cong (\mathbb{Z}_3 \times \mathbb{Z}_3) \rtimes \mathbb{Z}_3$ (order 27), generated by the clock and shift matrices $Q = \mathrm{diag}(1, \omega, \omega^2)$ and $P_{ij} = \delta_{i+1, j \bmod 3}$:
- **Universal 3-Generation Leptons & Quarks over $\mathbb{R}$:** For any self-conjugate spin/electroweak representation $\rho_W$ ($\iota(\rho_W) = \pm 1$), the Frobenius–Schur indicator of $\rho_W \boxtimes \mathbf{3}_f$ is:
  $$\iota(\rho_W \boxtimes \mathbf{3}_f) = \iota(\rho_W)\,\iota(\mathbf{3}_f) = (\pm 1) \times 0 = 0.$$
  Thus $(\rho_W \boxtimes \mathbf{3}_f) \oplus (\rho_W \boxtimes \bar{\mathbf{3}}_f)$ forms a **purely complex matrix block** $M_{3 d_W}(\mathbb{C})$ in $\mathbb{R}[G]$ containing **all 3 complex Dirac/Weyl generations** of leptons (and similarly for quarks in $\rho_W \boxtimes \mathbf{3}_c \boxtimes \mathbf{3}_f$).
- **Full $3 \times 3$ Yukawa Matrix Spanning & Non-Trivial $V_{\text{CKM}} \neq I_3$:** In $\Delta(27)_f$, the flavor bilinear decomposes into all **9 distinct 1-dimensional characters** of $\Delta(27)_f^{\text{ab}} \cong \mathbb{Z}_3 \times \mathbb{Z}_3$:
  $$\overline{\mathbf{3}}_f \otimes \mathbf{3}_f = \bigoplus_{(p, q) \in \mathbb{Z}_3 \times \mathbb{Z}_3} \chi_{p,q},$$
  corresponding to the 9 linearly independent Weyl–Heisenberg matrices $W_{p,q} = P^p Q^q$, which span the **entire complex matrix algebra $M_3(\mathbb{C})_{\text{flavor}}$**. Because $P Q = \omega Q P \neq Q P$, Higgs fields transforming in distinct $\chi_{p,q}$ channels generate **non-commuting** up and down mass matrices $[M_u M_u^\dagger, M_d M_d^\dagger] \neq 0$, naturally yielding realistic mass hierarchies and full CKM/PMNS mixing ($V_{\text{CKM}} \neq I_3, J_{\text{CP}} \neq 0$). $\blacksquare$

---

### Theorem 2.3 (Condition 3: Internal Chiral Gauge Completeness)
To cancel the $\mathrm{SU}(3)_c^3$, $\mathrm{SU}(2)_{W,L}^2 \mathrm{U}(1)_Y$, $\mathrm{U}(1)_Y^3$, and gravitational anomalies within each generation while predicting the bare unification Weinberg angle $\sin^2\theta_W(\Lambda) = 3/8$, the internal gauge group $G_{\text{gauge}}$ must contain:
- **Weak Left/Right Doublets:** $G_{W,L} \times G_{W,R}$ with quaternionic doublets $\mathbf{2}_{W,L}$ ($\iota = -1$) and $\mathbf{2}_{W,R}$ ($\iota = -1$) (minimal faithful $\mathrm{SU}(2)$ subgroup: $2T$ of order 24, or $Q_8$ of order 8), and
- **Color Triplets:** $G_c \subset \mathrm{SU}(3)_c$ with a faithful complex 3D triplet $\mathbf{3}_c$ ($\iota(\mathbf{3}_c) = 0$) and $\mathbb{Z}_3$ center $Z(G_c) \cong \mathbb{Z}_3$ whose 8 non-trivial central cosets admit a transitive action of $2T_W \cong \mathrm{SL}(2,\mathbb{F}_3)$ (minimal subgroup: $\Delta(27)_c$ of order 27 and exponent 3).

Under $G_{\text{gauge}} \times G_{\text{flavor}}$, a complete 3-generation Standard Model (Left–Right / Pati–Salam) fermion spectrum consists of:
$$\begin{aligned}
Q_L &\sim (\mathbf{2}_{\text{spin}, L},\; \mathbf{2}_{W,L},\; \mathbf{1}_{W,R},\; \mathbf{3}_c,\; \mathbf{3}_f), \qquad & Q_R &\sim (\mathbf{2}_{\text{spin}, R},\; \mathbf{1}_{W,L},\; \mathbf{2}_{W,R},\; \mathbf{3}_c,\; \mathbf{3}_f), \\
L_L &\sim (\mathbf{2}_{\text{spin}, L},\; \mathbf{2}_{W,L},\; \mathbf{1}_{W,R},\; \mathbf{1}_c,\; \mathbf{3}_f), \qquad & L_R &\sim (\mathbf{2}_{\text{spin}, R},\; \mathbf{1}_{W,L},\; \mathbf{2}_{W,R},\; \mathbf{1}_c,\; \mathbf{3}_f),
\end{aligned}$$
giving $N(\mathbf{3}_L) = N(\mathbf{3}_R) = 6$, exact cancellation of all gauge and gravitational anomalies, and $\sin^2\theta_W(\Lambda) = \frac{\mathrm{Tr}(T_{3L}^2)}{\mathrm{Tr}(T_{3L}^2) + \mathrm{Tr}(Y^2)} = \frac{3}{8}$. $\blacksquare$

---

### Theorem 2.4 (Condition 4: Abelianization Kirchhoff Law & Poset Washout Protection)
When $G = G_{\text{st}} \times G_{\text{gauge}} \times G_{\text{flavor}}$ includes factors with non-trivial abelianization (such as $2T^{\text{ab}} \cong \mathbb{Z}_3$ and $\Delta(27)^{\text{ab}} \cong \mathbb{Z}_3 \times \mathbb{Z}_3$), the total abelianization $G^{\text{ab}}$ is a non-trivial finite abelian 3-group $(\mathbb{Z}_3)^r$ ($r \ge 3$). By the Abelianization Kirchhoff Theorem ([`formalization.md`](formalization.md#L73-L94)) and Shared-Permutation Locking ([`falsification_report.md`](falsification_report.md#L119-L137)):
1. Every vertex of the poset $Z$ conserves $r$ independent additive $\mathbb{Z}_3$ charges ($\sum_{\text{in}} \pi_{\text{ab}}(u_i) \equiv \sum_{\text{out}} \pi_{\text{ab}}(v_j) \pmod 3$),
2. Even at high vertex degree $d \to \infty$, causal selectivity never drops below $S_\infty(G) = \log_2 \lvert G^{\text{ab}} \rvert = r \log_2 3 \ge 4.75\text{ bits}$, and
3. At manifold-like degrees $d \in \{3, 4\}$, shared permutations $(\sigma, \tau) \in S_d^2$ lock the perfect spacetime factor $G_{\text{st}} \cong 2I$ to the internal gauge/flavor factors, preventing causal washout of spacetime geometry. $\blacksquare$

---

## 3. Condition 5 (From Phase 3 & Phase 4 Audits): Hypercharge–Peccei–Quinn Separation

Cross-checking the Conjugate-Higgs Reality Trap ([`popperian_attacks_phase3.md`](popperian_attacks_phase3.md)) against the Electroweak Hypercharge / Weinberg Angle audit ([`constants_of_nature_phase4.md`](constants_of_nature_phase4.md)) reveals a **5th necessary condition**:

### Theorem 3.1 (Condition 5: Why Hypercharge $2T_{3R}$ and the Type-II 2HDM Peccei–Quinn Grading $\mathbb{Z}_3^{\text{PQ}}$ Cannot Share a Single $\mathbb{Z}_3$)
In $\Sigma(216 \times 3) \cong \Delta(27)_c \rtimes 2T_W$, the abelianization is a single $\mathbb{Z}_3$ ($2T_W^{\text{ab}} \cong \mathbb{Z}_3$), while bare $\Delta(27)_f$ possesses only **one** 3D flavor irrep $\Delta[9]$ of central flavor triality $+1$. Asking $\Sigma(216 \times 3)^{\text{ab}} \cong \mathbb{Z}_3$ to serve simultaneously as **continuous electroweak hypercharge $2T_{3R}$** and as the **discrete Type-II 2HDM Peccei–Quinn symmetry $\mathbb{Z}_3^{\text{PQ}}$** creates an inescapable contradiction (Attack 17):
1. **If $\Sigma(216 \times 3)^{\text{ab}} \cong \mathbb{Z}_3$ is used as continuous hypercharge / $2T_{3R}$** ($Q_L \sim \Sigma[18](0), u_R \sim \Sigma[11](+1), d_R \sim \Sigma[10](-1), L_L \sim \Sigma[5](0), e_R \sim \Sigma[2](-1), \nu_R \sim \Sigma[0](0)$), then $e_R$ is properly charged, $\nu_R$ is neutral, and the perturbative $[\mathrm{SU}(3)_c]^2 \mathrm{U}(1)_Y$ triangle anomaly vanishes in $\mathbb{R}$ ($2(0) - (+1) - (-1) = 0$). **However**, $H_u \sim \Sigma[4](\lambda^{+1})$ and $H_d \sim \Sigma[3](\lambda^{-1})$ become complex conjugates inside $M_2(\mathbb{C})$, forcing $\mathrm{Tr}(M_u M_u^\dagger) = \mathrm{Tr}(M_d M_d^\dagger)$ ($m_t = m_b$, Attack 8a) since the real octet $\Sigma[19](8_{+1})$ carries color ($\lvert\ker\Sigma[19]\rvert = 3 < 27$, Attack 8b).
2. **If $\Sigma(216 \times 3)^{\text{ab}} \cong \mathbb{Z}_3$ is instead used as the Type-II 2HDM Peccei–Quinn symmetry** (Attack 8c: $Q_L \sim \Sigma[18](0), u_R \sim \Sigma[7](0), d_R \sim \Sigma[10](-1), L_L \sim \Sigma[4](+1), \nu_R \sim \Sigma[1](+1), e_R \sim \Sigma[0](0)$), then $H_u = \Sigma[5] \in \mathbb{H}$ and $H_d = \Sigma[3] \in M_2(\mathbb{C})$ are separated into disjoint color-singlet blocks ($m_t \gg m_b$), **but**:
   - The right-handed electron $e_R \sim \Sigma[0]$ is placed in the **trivial representation** of $\Sigma(216 \times 3)$ ($\lvert\ker\Sigma[0]\rvert = 648$), making $e_R$ electrically neutral ($Q_{\text{em}}(e_R) = 0$) while the sterile neutrino $\nu_R \sim \Sigma[1](\lambda^{+1})$ is charged!
   - And $2q(Q_L) - q(u_R) - q(d_R) = +1 \neq 0$ leaves the continuous $[\mathrm{SU}(3)_c]^2 \mathrm{U}(1)_Y$ triangle diagram non-zero in $\mathbb{R}$ (since $3 \times 1 \equiv 0 \pmod 3$ cancels only *discrete* instanton anomalies, not perturbative continuous gauge anomalies).

**Minimal Resolution (Option D: $G = 2I \times \Sigma(216 \times 3) \times \Delta(27)_f \times \mathbb{Z}_3^{\text{PQ}}$, Order $6,298,560$):**
Extending the flavor sector by the discrete Peccei–Quinn factor $\mathbb{Z}_3^{\text{PQ}}$ splits $\mathbf{3}_f$ into three distinct 3D irreps $\mathbf{3}_f^{(0)}, \mathbf{3}_f^{(+1)}, \mathbf{3}_f^{(-1)}$ of central flavor charge $+1$, cleanly separating the two roles:
- $\Sigma(216 \times 3)$ carries **exact continuous gauge quantum numbers** ($Q_L \sim \Sigma[18](0), u_R \sim \Sigma[11](+1), d_R \sim \Sigma[10](-1), L_L \sim \Sigma[5](0), e_R \sim \Sigma[2](-1), \nu_R \sim \Sigma[0](0)$), and
- $\mathbb{Z}_3^{\text{PQ}}$ carries the **discrete Type-II 2HDM Peccei–Quinn charges** $q_{\text{PQ}}(Q_L, u_R, d_R, L_L, e_R, \nu_R) = (0, 0, -1, +1, 0, +1)$ (anomaly-free modulo 3 across the $N_{\text{gen}} = 3$ generations), placing $H_u = (\Sigma[4], \text{PQ}[0]) \in M_2(\mathbb{C})_{q_{\text{PQ}}=0}$, $H_d = H_e = (\Sigma[3], \text{PQ}[-1]) \in M_2(\mathbb{C})_{q_{\text{PQ}}=\mp 1}$, and $H_\nu = (\Sigma[5], \text{PQ}[0]) \in \mathbb{H}_{q_{\text{PQ}}=0}$ into **three completely disjoint, color-singlet Wedderburn blocks** of $\mathbb{R}[G]$. $\blacksquare$
