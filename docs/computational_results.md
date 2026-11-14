# Computational Exploration: Finite Group Candidates, Real Wedderburn Algebras, and Poset Conservation Laws

All results in this report are computed by the exact group-theoretic and spectral engine in [`src/`](../src/) via [`run_analysis.py`](../run_analysis.py).

---

## 1. Candidate Finite Groups and Structural Invariants

| Group $G$ | Order $\lvert G \rvert$ | Conjugacy Classes $k(G)$ | Commutator $\lvert [G,G] \rvert$ | Abelianization $\lvert G^{\text{ab}} \rvert$ | Perfect ($G=[G,G]$)? | Commuting Prob $\operatorname{cp}(G)$ |
| --- | --- | --- | --- | --- | --- | --- |
| **S_3** | 6 | 3 | 3 | 2 | No | 0.50000 (3/6) |
| **A_4** | 12 | 4 | 4 | 3 | No | 0.33333 (4/12) |
| **2T** | 24 | 7 | 8 | 3 | No | 0.29167 (7/24) |
| **Delta(27)** | 27 | 11 | 3 | 9 | No | 0.40741 (11/27) |
| **Delta(54)** | 54 | 10 | 27 | 2 | No | 0.18519 (10/54) |
| **Sigma(36x3)** | 108 | 14 | 27 | 4 | No | 0.12963 (14/108) |
| **2I** | 120 | 9 | 120 | 1 | Yes | 0.07500 (9/120) |
| **PSL(2, F_7)** | 168 | 6 | 168 | 1 | Yes | 0.03571 (6/168) |
| **Sigma(216)** | 216 | 10 | 72 | 3 | No | 0.04630 (10/216) |
| **Sigma(216x3)** | 648 | 24 | 216 | 3 | No | 0.03704 (24/648) |

---

## 2. Real Group Algebras $\mathbb{R}[G]$ and Connes' Standard Model Algebra $\mathcal{A}_F$

Using the Frobenius–Schur indicator $\iota(\chi) = \frac{1}{\lvert G \rvert} \sum_{g \in G} \chi(g^2) \in \{+1, 0, -1\}$, we decompose $\mathbb{R}[G]$ into simple matrix algebras over $\mathbb{R}$, $\mathbb{C}$, and $\mathbb{H}$ and test for the presence of Connes' internal Standard Model algebra $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$:

### 2.1 Single-Factor Groups

| Group $G$ | Irrep Degrees $d_\chi$ (with $\iota(\chi)$) | Real Wedderburn Decomposition $\mathbb{R}[G]$ | $\mathbb{H}$ Blocks ($d=2, \iota=-1$) | $M_3(\mathbb{C})$ Blocks ($d=3, \iota=0$) | $M_6(\mathbb{C})$ Blocks ($d=6, \iota=0$) |
| --- | --- | --- | --- | --- | --- |
| **S_3** | $1_{+1}$, $1_{+1}$, $2_{+1}$ | $R^2 (+) M_2(R)$ | 0 | 0 | 0 |
| **A_4** | $1_{+1}$, $1_{+0}$, $1_{+0}$, $3_{+1}$ | $R (+) M_3(R) (+) C$ | 0 | 0 | 0 |
| **2T** | $1_{+1}$, $1_{+0}$, $1_{+0}$, $2_{+0}$, $2_{+0}$, $2_{-1}$, $3_{+1}$ | $R (+) M_3(R) (+) C (+) M_2(C) (+) H$ | 1 | 0 | 0 |
| **Delta(27)** | $1_{+1}$, $1_{+0}$, $1_{+0}$, $1_{+0}$, $1_{+0}$, $1_{+0}$, $1_{+0}$, $1_{+0}$, $1_{+0}$, $3_{+0}$, $3_{+0}$ | $R (+) C^4 (+) M_3(C)$ | 0 | 1 | 0 |
| **Delta(54)** | $1_{+1}$, $1_{+1}$, $2_{+1}$, $2_{+1}$, $2_{+1}$, $2_{+1}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$ | $R^2 (+) M_2(R)^4 (+) M_3(C)^2$ | 0 | 2 | 0 |
| **Sigma(36x3)** | $1_{+1}$, $1_{+1}$, $1_{+0}$, $1_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $4_{+1}$, $4_{+1}$ | $R^2 (+) M_4(R)^2 (+) C (+) M_3(C)^4$ | 0 | 4 | 0 |
| **2I** | $1_{+1}$, $2_{-1}$, $2_{-1}$, $3_{+1}$, $3_{+1}$, $4_{+1}$, $4_{-1}$, $5_{+1}$, $6_{-1}$ | $R (+) M_3(R)^2 (+) M_4(R) (+) M_5(R) (+) H^2 (+) M_2(H) (+) M_3(H)$ | 2 | 0 | 0 |
| **PSL(2, F_7)** | $1_{+1}$, $3_{+0}$, $3_{+0}$, $6_{+1}$, $7_{+1}$, $8_{+1}$ | $R (+) M_6(R) (+) M_7(R) (+) M_8(R) (+) M_3(C)$ | 0 | 1 | 0 |
| **Sigma(216)** | $1_{+1}$, $1_{+0}$, $1_{+0}$, $2_{+0}$, $2_{+0}$, $2_{-1}$, $3_{+1}$, $8_{+1}$, $8_{+0}$, $8_{+0}$ | $R (+) M_3(R) (+) M_8(R) (+) C (+) M_2(C) (+) M_8(C) (+) H$ | 1 | 0 | 0 |
| **Sigma(216x3)** | $1_{+1}$, $1_{+0}$, $1_{+0}$, $2_{+0}$, $2_{+0}$, $2_{-1}$, $3_{+1}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $3_{+0}$, $6_{+0}$, $6_{+0}$, $6_{+0}$, $6_{+0}$, $6_{+0}$, $6_{+0}$, $8_{+1}$, $8_{+0}$, $8_{+0}$, $9_{+0}$, $9_{+0}$ | $R (+) M_3(R) (+) M_8(R) (+) C (+) M_2(C) (+) M_3(C)^3 (+) M_6(C)^3 (+) M_8(C) (+) M_9(C) (+) H$ | 1 | 3 | 3 |

### 2.2 Three Major Algebraic Discoveries from Exact Computation

1. **Correction to the $2I$ Irrep Spectrum & Why $2I$ Alone Cannot Produce $M_3(\mathbb{C})$:**
   - $2I \cong \operatorname{SL}(2, \mathbb{F}_5)$ (order 120) has **9** irreducible representations of dimensions $1_{+1}, 2_{-1}, 2_{-1}, 3_{+1}, 3_{+1}, 4_{+1}, 4_{-1}, 5_{+1}, 6_{-1}$ ($\sum d_\chi^2 = 120$; note that Exchange 16/18 of the Gemini chat omitted one of the two 4D irreps).
   - Every irrep of $2I$ has Frobenius–Schur indicator $\iota = +1$ (real, factoring through $\operatorname{SO}(3) \supset A_5$) or $\iota = -1$ (quaternionic, faithful in $\operatorname{SU}(2)$):
     $$\mathbb{R}[2I] \cong \mathbb{R} \oplus M_3(\mathbb{R})^2 \oplus M_4(\mathbb{R}) \oplus M_5(\mathbb{R}) \oplus \mathbb{H}^2 \oplus M_2(\mathbb{H}) \oplus M_3(\mathbb{H}).$$
   - Because $2I$ has **zero** complex irreps ($\iota = 0$), its 3D irreps generate $M_3(\mathbb{R})$ rather than $M_3(\mathbb{C})$. Thus $2I$ natively supplies the **quaternionic electroweak doublets ($\mathbb{H}^2$)**, but **must** be paired with a complex-triality group ($\iota(\chi_3) = 0$) to yield $\operatorname{SU}(3)_c$ color.
2. **Projective $\Sigma(216) \subset \operatorname{PSU}(3)$ vs. Faithful $\Sigma(216 \times 3) \subset \operatorname{SU}(3)$:**
   - The centerless Hessian group $\Sigma(216) \cong (\mathbb{Z}_3 \times \mathbb{Z}_3) \rtimes \operatorname{SL}(2, \mathbb{F}_3)$ of order 216 is a subgroup of $\operatorname{PSU}(3) = \operatorname{SU}(3)/\mathbb{Z}_3$. Because fundamental 3D triplets transform faithfully under the center $Z(\operatorname{SU}(3)) \cong \mathbb{Z}_3$, **$\mathbb{R}[\Sigma(216)]$ contains zero $M_3(\mathbb{C})$ blocks** (its 3D irrep is real, $M_3(\mathbb{R})$, and its complex irreps have degrees $1, 2, 8$).
   - By contrast, its triple cover **$\Sigma(216 \times 3) \subset \operatorname{SU}(3)$ of order 648** has Real Wedderburn decomposition:
     $$\mathbb{R}[\Sigma(216 \times 3)] \cong \mathbb{R} \oplus M_3(\mathbb{R}) \oplus M_8(\mathbb{R}) \oplus \mathbb{C} \oplus M_2(\mathbb{C}) \oplus M_3(\mathbb{C})^3 \oplus M_6(\mathbb{C})^3 \oplus M_8(\mathbb{C}) \oplus M_9(\mathbb{C}) \oplus \mathbb{H}.$$
   - Look at what sits inside this single group of order 648: Connes' $\mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ directly, **plus exactly 3 generations of $M_3(\mathbb{C})$ color triplets and 3 generations of $M_6(\mathbb{C})$ electroweak-color $(2, 3)$ quark multiplets!**
3. **Comparison of $\operatorname{SU}(3)$ Triality Series ($\Delta(27) \subset \Delta(54) \subset \Sigma(36 \times 3) \subset \Sigma(216 \times 3)$) and $\operatorname{PSL}(2, \mathbb{F}_7)$:**
   - $\operatorname{PSL}(2, \mathbb{F}_7)$ (order 168) and $\Delta(27)$ (order 27) each provide **1** copy of $M_3(\mathbb{C})$.
   - $\Delta(54)$ provides **2** copies of $M_3(\mathbb{C})$.
   - $\Sigma(36 \times 3)$ (order 108) provides **4** copies of $M_3(\mathbb{C})$.
   - $\Sigma(216 \times 3)$ (order 648) provides **3** copies of $M_3(\mathbb{C})$—matching the observed **3 fermion generations ($N_g = 3$)**.

### 2.3 Composite Product Candidates ($G_1 \times G_2$)

Using $\iota(\chi_1 \boxtimes \chi_2) = \iota(\chi_1)\iota(\chi_2)$, the product candidates decompose as follows:

| Composite Group $G$ | Order $\lvert G \rvert$ | $1\text{D}$ Scalar Blocks ($\mathbb{R}$ or $\mathbb{C}$) | $\mathbb{H}$ Weak Blocks ($d=2, \iota=-1$) | $M_3(\mathbb{C})$ Color Blocks ($d=3, \iota=0$) | $M_6(\mathbb{C})$ Quark Bi-Fundamental Blocks |
| --- | --- | --- | --- | --- | --- |
| **2I x Delta(27)** | 3,240 | 5 | 2 | 9 | 6 |
| **2I x Delta(54)** | 6,480 | 2 | 4 | 2 | 4 |
| **2I x Sigma(36x3)** | 12,960 | 3 | 4 | 6 | 9 |
| **2I x PSL(2, F_7)** | 20,160 | 1 | 2 | 1 | 2 |
| **2I x Sigma(216)** | 25,920 | 2 | 3 | 2 | 3 |
| **2I x Sigma(216x3)** | 77,760 | 2 | 3 | 5 | 12 |

---

## 3. Cayley Graph Spectra and Bare Higgs / Gauge Coupling Ratios

For each symmetric generating set $S = S^{-1} \subset G \setminus \{e\}$, we compute the exact 3-cycle count $W_3(S)$, 4-cycle (plaquette) count $W_4(S)$, spectral quartic ratios $\kappa(A_S) = \lvert G \rvert \operatorname{Tr}(A_S^4) / (\operatorname{Tr}(A_S^2))^2 = W_4(S)/\lvert S \rvert^2$ and $\kappa(\Delta_S)$, and the resulting tree-level bare ratio $m_H(\Lambda)/m_W(\Lambda) = 2\sqrt{\kappa}$:

| Group $G$ | Generating Set $S$ | Degree $\lvert S \rvert$ | Triangles $W_3(S)$ | Plaquettes $W_4(S)$ | $\kappa(A_S)$ | $\kappa(\Delta_S)$ | $2\sqrt{\kappa(\Delta_S)}$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **2I** | 600-cell edges (C_12, ord 10) | 12 | 60 | 636 | 4.4167 | 1.1859 | 2.1780 |
| **2I** | Icosahedral reflections (C_30, ord 4) | 30 | 120 | 8520 | 9.4667 | 1.1170 | 2.1138 |
| **2I** | Dodecahedral faces (C_20, ord 6) | 20 | 120 | 2220 | 5.5500 | 1.1373 | 2.1329 |
| **PSL(2, F_7)** | Involutions (C_21, ord 2) | 21 | 84 | 1953 | 4.4286 | 1.1476 | 2.1425 |
| **PSL(2, F_7)** | Octonionic 4-cycles (C_42, ord 4) | 42 | 672 | 23016 | 13.0476 | 1.0628 | 2.0618 |
| **PSL(2, F_7)** | Fano 7-cycles (C_24 + C_24*, deg 48) | 48 | 576 | 33408 | 14.5000 | 1.0656 | 2.0646 |
| **Delta(27)** | Minimal Heisenberg gens (deg 4) | 4 | 4 | 28 | 1.7500 | 1.5100 | 2.4576 |
| **Delta(27)** | All non-central classes (deg 24) | 24 | 504 | 12312 | 21.3750 | 1.0518 | 2.0511 |
| **Sigma(216)** | Involutions (C_9, ord 2) | 9 | 0 | 729 | 9.0000 | 1.4400 | 2.4000 |
| **Sigma(216x3)** | Self-inverse 3-cycles (C_24, ord 3) | 24 | 504 | 12312 | 21.3750 | 1.0518 | 2.0511 |
| **Sigma(216x3)** | Order-4 class (C_54, ord 4) | 54 | 1080 | 49896 | 17.1111 | 1.0503 | 2.0497 |

### 3.1 Exact Spectrum of the 600-Cell Cayley Graph ($2I$, $\lvert S \rvert = 12$)

The 12 nearest neighbors of the identity in the 600-cell form a single conjugacy class of order-10 elements in $2I$. Via Schur's lemma ($\lambda_\chi = \frac{12 \chi(C_{12})}{d_\chi}$ with multiplicity $d_\chi^2$), the exact adjacency spectrum on $\ell^2(2I)$ has **9** distinct eigenvalues (correcting the 6-value list in the Gemini chat):

$$\operatorname{Spec}(A_{600}) = \{ 12.0000^{(1)}, 9.7082^{(4)}, 6.4721^{(9)}, 3.0000^{(16)}, 0.0000^{(25)}, -2.0000^{(36)}, -2.4721^{(9)}, -3.0000^{(16)}, -3.7082^{(4)} \}$$

in exact agreement with $\{12^{(1)}, (6\phi)^{(4)}, (4\phi)^{(9)}, 3^{(16)}, 0^{(25)}, (-2)^{(36)}, (-3)^{(16)}, (-4\phi^{-1})^{(9)}, (-6\phi^{-1})^{(4)}\}$.

---

## 4. Ordered Product Conservation Law ($\Pi(\mathbf{g}_{\text{in}}) \cap \Pi(\mathbf{g}_{\text{out}}) \neq \emptyset$)

### 4.1 Exact Vertex Compatibility Rates and Closed-Form Verification

| Group $G$ | $\lvert G \rvert$ | $\mathbb{E}[\lvert\Pi(g_1, g_2)\rvert]$ | Exact $p_{1,2}(G)$ | Formula $\frac{2\lvert G \rvert - k(G)}{\lvert G \rvert^2}$ | Exact $p_{2,2}(G)$ | $\mathbb{E}[\lvert\Pi(g_1,g_2,g_3)\rvert]$ | Est. $p_{3,3}(G)$ | $\mathbb{E}[\lvert\Pi(g_1..g_4)\rvert]$ | Est. $p_{4,4}(G)$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **S_3** | 6 | 1.5000 | 0.250000 | 0.250000 | 0.319444 | 2.16 | 0.4453 | 2.57 | 0.4887 |
| **A_4** | 12 | 1.6667 | 0.138889 | 0.138889 | 0.204861 | 2.95 | 0.3148 | 3.61 | 0.3303 |
| **2T** | 24 | 1.7083 | 0.071181 | 0.071181 | 0.110243 | 3.37 | 0.2548 | 5.56 | 0.3171 |
| **Delta(27)** | 27 | 1.5926 | 0.058985 | 0.058985 | 0.080933 | 2.46 | 0.1049 | 2.83 | 0.1116 |
| **Delta(54)** | 54 | 1.8148 | 0.033608 | 0.033608 | 0.057613 | 4.22 | 0.2304 | 10.44 | 0.4399 |
| **Sigma(36x3)** | 108 | 1.8704 | 0.017318 | 0.017318 | 0.030907 | 4.63 | 0.1353 | 11.68 | 0.2414 |
| **2I** | 120 | 1.9250 | 0.016042 | 0.016042 | 0.030097 | 5.18 | 0.1864 | 16.29 | 0.7948 |
| **PSL(2, F_7)** | 168 | 1.9643 | 0.011692 | 0.011692 | 0.022659 | 5.60 | 0.1654 | 19.36 | 0.8464 |
| **Sigma(216)** | 216 | 1.9537 | 0.009045 | 0.009045 | 0.017333 | 5.44 | 0.1101 | 17.32 | 0.3119 |
| **Sigma(216x3)** | 648 | 1.9630 | 0.003029 | 0.003029 | 0.005859 | 5.58 | 0.0442 | 19.14 | 0.2554 |

### 4.2 Causal Selectivity vs. Vertex Degree: Why Bounded-Degree Posets Beat Kleitman–Rothschild Posets

For each vertex degree $d = \lvert\operatorname{Pred}_Z(z)\rvert = \lvert\operatorname{Succ}_Z(z)\rvert$, we measure the commutator-coset fill fraction $\mathbb{E}[\lvert\Pi(\mathbf{g})\rvert]/\lvert[G,G]\rvert$, the $(d \to d)$ compatibility probability $p_{d,d}(G)$, and the **causal selectivity** $S_d(G) = -\log_2 p_{d,d}(G)$ (in bits of constraint information per vertex):

| Group $G$ | $d=1$ ($S_1$ bits) | $d=2$ ($S_2$ bits) | $d=3$ ($S_3$ bits) | $d=4$ ($S_4$ bits) | $d=5$ ($S_5$ bits) | $d=6$ ($S_6$ bits) |
| --- | --- | --- | --- | --- | --- | --- |
| **S_3** | 2.58b ($p=0.167$) | 1.65b ($p=0.319$) | 1.17b ($p=0.445$) | 1.03b ($p=0.489$) | 1.01b ($p=0.497$) | 1.00b ($p=0.500$) |
| **A_4** | 3.58b ($p=0.083$) | 2.29b ($p=0.205$) | 1.67b ($p=0.315$) | 1.60b ($p=0.330$) | 1.59b ($p=0.333$) | 1.59b ($p=0.333$) |
| **2T** | 4.58b ($p=0.042$) | 3.18b ($p=0.110$) | 1.97b ($p=0.255$) | 1.66b ($p=0.317$) | 1.59b ($p=0.332$) | 1.58b ($p=0.334$) |
| **Delta(27)** | 4.75b ($p=0.037$) | 3.63b ($p=0.081$) | 3.25b ($p=0.105$) | 3.16b ($p=0.112$) | 3.17b ($p=0.111$) | 3.16b ($p=0.112$) |
| **Delta(54)** | 5.75b ($p=0.019$) | 4.12b ($p=0.058$) | 2.12b ($p=0.230$) | 1.18b ($p=0.440$) | 1.04b ($p=0.488$) | 1.01b ($p=0.497$) |
| **Sigma(36x3)** | 6.75b ($p=0.009$) | 5.02b ($p=0.031$) | 2.89b ($p=0.135$) | 2.05b ($p=0.241$) | 2.00b ($p=0.250$) | 1.99b ($p=0.251$) |
| **2I** | 6.91b ($p=0.008$) | 5.05b ($p=0.030$) | 2.42b ($p=0.186$) | 0.33b ($p=0.795$) | 0.01b ($p=0.994$) | 0.00b ($p=1.000$) |
| **PSL(2, F_7)** | 7.39b ($p=0.006$) | 5.46b ($p=0.023$) | 2.60b ($p=0.165$) | 0.24b ($p=0.846$) | 0.00b ($p=0.999$) | -0.00b ($p=1.000$) |
| **Sigma(216)** | 7.75b ($p=0.005$) | 5.85b ($p=0.017$) | 3.18b ($p=0.110$) | 1.68b ($p=0.312$) | 1.59b ($p=0.332$) | 1.59b ($p=0.333$) |
| **Sigma(216x3)** | 9.34b ($p=0.002$) | 7.42b ($p=0.006$) | 4.50b ($p=0.044$) | 1.97b ($p=0.255$) | 1.59b ($p=0.332$) | 1.59b ($p=0.333$) |

### 4.3 Physical Interpretation of the Degree Sweep

1. **Degenerate Chains ($d=1$):** At $d=1$, $\lvert\Pi(g)\rvert = 1$, so no non-abelian commutator $[g_i, g_j]$ can ever be generated; the state is frozen along 1D chains.
2. **The Manifold-Like Window ($d \in \{2, 3, 4\}$):** For candidate groups like $2I$, $\operatorname{PSL}(2, \mathbb{F}_7)$, $\Sigma(216)$, and $\Sigma(216 \times 3)$, degrees $d \in \{2, 3, 4\}$ simultaneously generate rich non-abelian permutation sets ($\mathbb{E}[\lvert\Pi\rvert] \in [1.9, 19.3]$) **while maintaining high causal selectivity** ($1.7$ to $7.4$ bits of constraint per vertex).
3. **Kleitman–Rothschild Washout ($d \ge 5\text{--}6$):** In pathological 3-layer Kleitman–Rothschild posets, interior vertices have macroscopic degree $d \sim N/4$. As shown in the table above, once $d \ge 5\text{--}6$, $\Pi(g_1, \dots, g_d)$ rapidly fills the entire commutator coset $g_1 \cdots g_d [G,G]$. For perfect groups ($2I$ and $\operatorname{PSL}(2, \mathbb{F}_7)$ where $[G,G]=G$), $p_{d,d}(G) \to 1$ and the causal selectivity $-\log_2 p_{d,d}(G) \to 0.00\text{ bits}$—meaning high-degree Kleitman–Rothschild hubs suffer **complete causal washout** (zero mutual information between past and future).
